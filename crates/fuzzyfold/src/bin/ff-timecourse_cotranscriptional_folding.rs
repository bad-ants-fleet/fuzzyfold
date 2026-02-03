use clap::Parser;
use anyhow::Result;
use colored::*;
use serde_json::to_string_pretty;
use std::fs;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;
use rayon::prelude::*;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use rand::rng;

use ff_structure::PairTable;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_kinetics::Metropolis;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::timeline_pairlists::Timeline;
use ff_kinetics::timeline_plotting_pairlists::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistryPL;
//use ff_kinetics::timeline::Timeline;
//use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
//use ff_kinetics::MacrostateRegistry;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimelineParameters;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    /// Backup/Store timeline in this file.
    #[arg(long, value_name = "FILE")]
    timeline: Option<PathBuf>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0);

    let (header, sequence, _structure) = read_fasta_like_input(&cli.input)?;
    
    let name = if let Some(h) = header {
        println!("{}", h.yellow());
        h.strip_prefix('>')
            .and_then(|s| s.split_whitespace().next())
            .unwrap_or("anonymous")
            .to_string()
    } else {
        println!("{}", ">anonymous".yellow());
        "anonymous".to_string()
    };
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    // --Output times--

    let sequence_len = sequence.len();

    let t_lin = cli.simulation.t_lin; //timepoints per nucleotide
    let t_ext = cli.simulation.t_ext; //extension time (simulation time per nucleotide)
    let t_log = cli.simulation.t_log; //timepoints for posttranscriptional folding 
    let t_end = cli.simulation.t_end; //simulation stop time (Here: posttranscriptional stop time) (Only posttranscriptional time or total time??)
    
    let mut times = vec![0.0];

    // Linear segments: append `t_lin` evenly spaced points
    
    let mut len = 1;

    while len < sequence_len {
        let start = *times.last().unwrap();
        let step = t_ext / t_lin as f64;
        for i in 1..= t_lin {
            times.push(start + i as f64 * step);
        }
        len += 1;
    }

    // Logarithmic tail
    let start = *times.last().unwrap();
    let log_start = start.ln();
    let log_end = (start + t_end).ln();
    for i in 1..t_log {
        let frac = i as f64 / t_log as f64;
        let value = (log_start + frac * (log_end - log_start)).exp();
        times.push(value);
    }
    times.push(start + t_end);


    let mut registry = MacrostateRegistryPL::from((&sequence, &emodel));
    //let mut registry = MacrostateRegistry::from((&sequence, &emodel));
    let _ = registry.insert_files(&cli.macrostates);

    println!("Macrostates:\n{}", registry.iter()
        .map(|(_, m)| format!(" - {} {:6.2}", m.name(), m.ensemble_energy().unwrap_or(0.0)))
        .collect::<Vec<_>>().join("\n"));

    let shared_registry = Arc::new(registry);

    // If timeline.json exists, reload instead of starting empty
    let mut master = if let Some(path) = &cli.timeline {
        if Path::new(path).exists() {
            println!("Loading existing timeline from: {}", path.display());
            Timeline::from_file(path, &times, Arc::clone(&shared_registry))?
        } else {
            println!("A new timeline file will be created: {}", 
                path.display());
            Timeline::new(&times, Arc::clone(&shared_registry))
        }
    } else {
        Timeline::new(&times, Arc::clone(&shared_registry))
    };

    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    let timelines: Vec<_> = (0..cli.num_sims)
        .into_par_iter()
        .map_init(
            || pb.clone(), // each thread gets a clone
            |pb, _| {
                let registry = Arc::clone(&shared_registry);
                let mut timeline = Timeline::new(&times, registry);

                let mut t = 0.0;

                // start with len=1
                let mut current_len = 1;
                let mut current_struct = DotBracketVec(vec![DotBracket::Unpaired]); //initialize current structure with length 1
                let mut current_seq = &sequence[..current_len];

                let mut t_idx = 0;


                while current_len < sequence.len() { //repeat until end of transcription

                    current_seq = &sequence[..current_len];

                    let pairings = PairTable::try_from(&current_struct).unwrap(); //build PairTable
                    let loops = LoopStructure::try_from((&current_seq[..], &pairings, &emodel)).unwrap(); //build loop structure 

                    // --- Check if code panics, because no possible transitions for current seequence, skip simulation in this case ---
                    
                    let add_pair = loops // can a base pair be added?
                        .get_add_neighbors_per_loop()
                        .iter()
                        .any(|(_, add_neighbors) | !add_neighbors.is_empty());

                    let del_pair = !loops.get_del_neighbors().is_empty(); //can a base pair be deleted?

                    if !add_pair && !del_pair { // no possible transitions => continue without simulation and extend structure and sequence
                        for i in t_idx..(t_idx + t_lin){
                            timeline.assign_structure(i, &current_struct);
                        }
                        t_idx += t_lin;
                        t = t + t_ext;
                        current_len += 1;
                        current_struct.0.push(DotBracket::Unpaired);
                        continue;
                    }

                    let mut simulator = LoopStructureSSA::from((loops, &rmodel)); //build simulator from loop structure and rate model
                    
                    let mut final_struct = current_struct.clone();
                  
                    let idx = t_idx + t_lin;

                    simulator.simulate(
                        &mut rng(), //random number 
                        t_ext,
                        |t_sim, tinc, _flux, ls| {
                            while t_idx < idx + t_lin + 1 && t + t_sim + tinc >= times[t_idx] {
                                let structure = DotBracketVec::from(ls);
                                timeline.assign_structure(t_idx, &structure);
                                t_idx += 1;
                            }
                            final_struct = DotBracketVec::from(ls);
                            true
                        },
                    );

                    current_struct = final_struct.clone();
                    t = t + t_ext; //update t

                    // ---Append sequence with next nucleotide---
                    current_len += 1; 
                    current_seq = &sequence[..current_len];
                    current_struct.0.push(DotBracket::Unpaired);

                }


                // ---Postranscriptional folding---


                let pairings = PairTable::try_from(&current_struct).unwrap();//build PairTable
                let loops = LoopStructure::try_from((&current_seq[..], &pairings, &emodel)).unwrap(); //build loop structure

                let mut simulator = LoopStructureSSA::from((loops, &rmodel)); //build simulator from loop structure and rate model
                
                let cofolding_time = (sequence.len() -1) as f64 * t_ext;

                simulator.simulate(
                    &mut rng(),
                    cofolding_time + t_end,
                    |t_sim, tinc, _, ls| {
                        while t_idx < times.len() && t + t_sim + tinc >= times[t_idx] {
                            let structure = DotBracketVec::from(ls);
                            timeline.assign_structure(t_idx, &structure);
                            t_idx += 1;
                        }
                        true
                    },
                );
                pb.inc(1);
                timeline
            }
        ).collect();

    pb.finish_with_message("All simulations complete!");

    // Master timeline
    for timeline in timelines {
        master.merge(timeline);
    }

    let t_final = master.points.last().unwrap().time;

    let cofolding_time = (sequence.len() -1) as f64 * cli.simulation.t_ext; 

    println!("Final Timeline:\n{}", master);
    plot_occupancy_over_time(&master, &format!("ff_{}.svg", name), cofolding_time, t_final);

    if let Some(path) = cli.timeline {
        let serial = master.to_serializable();
        let json = to_string_pretty(&serial)?;
        fs::write(path, json)?;
    }

    for i in times {
        println!("{}", i);
    }

    Ok(())
}
