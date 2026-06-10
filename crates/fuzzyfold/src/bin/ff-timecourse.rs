//! Stochastic folding simulator 
//! 
//! This binary performs stochastic folding simulations both at full sequence length and 
//! cotranscriptional, running multiple simulations in parallel and merging them into a master timeline. 
//! This timeline is used to produce an occupancy plot. The plot shows the occupancy of predefined 
//! macrostates in average over time. 
//! 
//! Mode: 
//! Full length simulation: give structure of full sequence length 
//! Cotranscriptonal simulation: give no structure or structure shorter than sequence length
//! 
//! --- Full length simulation --- 
//! 
//! Parameters:
//!  
//! - t-end: total simulation time 
//! - t-sep: last timepoint on linear scale (timepoint where linear scale in occupany plot switches to logarithmic scale and seperator is placed)
//! - t-lin: timepoints on linear time scale [0...t-sep]
//! - t-log: time points on logarithmic timescale [t-sep...t-log]
//! - num-sims: number of simulations performed
//! - k0: kinetic rate constant 
//! 
//! Input
//! - user has to give structure and t-sep
//! - optional: t-lin and t-log (default: 100)
//! - optional: t-end (default: 1.0)
//! - optional: num-sim (default: 1)
//! 
//! --- Cotranscriptional simulation ---
//! 
//! Parameters
//!  
//! - t-ext: extension time (simulation time per transcript length)
//! - t-end: time of postranscriptional simulation
//! - t-sep: last timepoint on linear scale (timepoint where linear scale in occupany plot switches to logarithmic scale and seperator is placed)
//! - t-lin: recorded time steps per transcript length [without t-sep]/ timepoints recorded on linear timescale [with t-sep]
//! - t-log: time points recorded for the posttranscriptional simulation on a logarithmic timescale [without t-sep]/ timepoints recorded on logarithmic timescale [with t-sep]
//! - num-sims: number of simulations performed
//! - k0: kinetic rate constant 
//! 
//! Input & mode 
//! - Cotranscriptional simulation 
//! -> user has to give either no structure or a structure, that is shorter than full length (start structure) and t-ext (extension time)
//! - optional: t-sep 
//! -> no t-sep: linear scale ends at end of transcription (t-lin: timepoints recorded per extension step; t-log: timepoints recorded during posttranscriptional folding)
//! -> t-sep: linear timescale ends at t-sep (t-lin: timepoints recorded on linear timescale; t-log: time points recorded on logarithmic timescale)
//! - optional: t-lin and t-log (default: 100)
//! - optional: t-end (default: 1.0)
//! - optional: num-sim (default: 1)
//! 
//! --- General ---
//! 
//! Output: 
//! - user has to give output file
//! - output formats: 
//! -> svg (occupancy plot)
//! -> nxy (occupancies)
//! -> tln (timeline)
//! 
//! Load timeline:
//! The 'tln' file can be used to accumulate results incrementially. Existing timlines can be reloaded and when a simulation is
//! run with it, it gets updated. (Parameters have to be set to the same values!)
//! 



use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;

use plotters::prelude::LogScalable;
use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use serde_json::to_string_pretty;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;
use ff_kinetics::timeline_pl::Timeline;
use ff_kinetics::timeline_plotting_pl::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistryPL;

use fuzzyfold::input_parsers::read_cotr_input;
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

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_cotr_input(&cli.input, is_rna)?;
    let pairings = PairTable::try_from(&structure)?;

    // --- Check input ---
    if structure.len() < sequence.len() && cli.simulation.t_ext.is_none() {
        panic!("Error:'t-ext' (extension time) missing for cotranscriptional simulation!"); 
    }

    if structure.len() == sequence.len() && cli.simulation.t_sep.is_none() {
        panic!("Error: 't-sep' (last timepoint of the linear scale) missing for full length simulation!")
    }

    if structure.len() == sequence.len() && cli.simulation.t_ext.is_some() {
        panic!("Error: 't-ext' (extension time) given in combination with full length structure. For cotranscriptional simulation give either 
        shorter start structure or no structure!");
    }

    if cli.simulation.t_ext.is_some() && cli.simulation.t_sep.is_some() {
        if cli.simulation.t_sep.unwrap() >= (cli.simulation.t_ext.unwrap() * (sequence.len() - structure.len()).as_f64() + cli.simulation.t_end) {
            panic!("Error: 't_sep' must be smaller than the total simulation time!");
        } 
    }

    let n = header.clone();
    let name = n
        .as_deref()
        .map(|h| h.trim_start_matches('>'))
        .and_then(|h| h.split_whitespace().next())
        .unwrap_or("unnamed_sequence");

    if let Some(h) = header {
        println!("{}", h.yellow());
    }

    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times(&sequence, Some(&structure)); // build times vector for output times 

    let mut macrostates = MacrostateRegistryPL::from((sequence.clone(), emodel.clone()));
    macrostates.insert_files(&cli.macrostates)?;
    // Verbose Output
    println!("{:>4} {:<10} {} {:>5} {:>8}",
        "ID",
        "Macrostate".cyan(), format!("{}", sequence).yellow(), "Size", "Energy");
    for (id, m) in macrostates.iter() {
        if m.name() == "Unassigned" {
            println!("{:4} {:<10}", 0, m.name());
            continue
        }
        println!("{:4} {:<10} {:<} {:>5} {:>8.2}",
            id, 
            m.name(),
            m.get_lowest_microstate().unwrap(),
            m.len(),
            m.ensemble_energy().unwrap());
    }
    let shared_macrostates = Arc::new(macrostates);


    // Output paths 
    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");
    let nxy_path = cli.output.with_extension("nxy");

    // If timeline.json exists, reload instead of starting empty
    let mut master = 
        if Path::new(&tln_path).exists() {
            println!("Loading existing timeline from: {}", tln_path.display());
            Timeline::from_file(&tln_path, &times, Arc::clone(&shared_macrostates))?
        } else {
            println!("A new timeline file will be created: {}", tln_path.display());
            Timeline::new(&times, Arc::clone(&shared_macrostates))
        };


    if cli.simulation.t_ext.is_none() { // full length simulation 
        let timelines: Vec<_> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
        };

        for timeline in timelines {
            master.merge(timeline);
        }
        
        plot_occupancy_over_time(&master, svg_path.clone(), &format!("Timecourse {}", name),cli.simulation.t_sep.unwrap(), cli.simulation.t_end);


    } else { // cotranscriptional folding simulation

        let mut sim_times = vec![cli.simulation.t_ext.unwrap(); sequence.len() - structure.len()];
        sim_times.push(cli.simulation.t_end);

        let timelines: Vec<_> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_cotimecourse(moves, rmodel, &sim_times, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
        };

        for timeline in timelines {
            master.merge(timeline);
        }

        // cotranscriptional simulation time
        let co_time = cli.simulation.t_ext.unwrap() * (sequence.len() - structure.len()).as_f64(); 
        // total simulation time
        let post_time = co_time + cli.simulation.t_end;

        if cli.simulation.t_sep.is_none() { // No t-sep given
            
            plot_occupancy_over_time(&master, svg_path.clone(), &format!("Cotimecourse {}", name),co_time, post_time);

        } else { // t-sep given 

            plot_occupancy_over_time(&master, svg_path.clone(), &format!("Cotimecourse {}", name), cli.simulation.t_sep.unwrap(), post_time);  
        }

    }


    println!("{}", "Finished simulations!".red());

    // save / print / plot.
    let mut writer = BufWriter::new(File::create(nxy_path.clone())?);
    write!(writer, "{}", master)?;
    println!("Wrote nxy file: {}", format!("{}",nxy_path.display()).green());
    
    println!("Plotted svg file: {}", svg_path.display());
    let serial = master.to_serializable();
    let json = to_string_pretty(&serial).unwrap();
    fs::write(tln_path.clone(), json).unwrap();
    println!("Wrote tln file: {}", tln_path.display());

    Ok(())
}


/// full length folding simulation 
fn run_timecourse<W, K, E>(
    moves: W,
    rmodel: K,
    t_end: f64,
    num_sims: u64,
    registry: Arc<MacrostateRegistryPL<E>>,
    times: &[f64],
) -> impl ParallelIterator<Item = Timeline<E>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
    E: EnergyModel + Send + Sync,
    SSA<W, K>: From<(W, K)>,
{
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _| {
                let registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(times, registry);

                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let mut t_idx = 0;
                simulator.simulate(
                    &mut rng(),
                    t_end,
                    |t, tinc, _, w| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            let structure = w.current_structure();
                            timeline.assign_structure(t_idx, &structure);
                            t_idx += 1;
                        }
                        true
                    },
                );

                pb.inc(1);
                timeline
            },
        )
}


/// Cotranscriptional folding simulation 
fn run_cotimecourse<W, K, E>(
    moves: W,
    rmodel: K,
    sim_times: &[f64],
    num_sims: u64,
    registry: Arc<MacrostateRegistryPL<E>>,
    times: &[f64],
) -> impl ParallelIterator<Item = Timeline<E>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
    E: EnergyModel + Send + Sync,
    SSA<W, K>: From<(W, K)>,
{
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    (0..num_sims)
    .into_par_iter()
    .map_init(
        move || pb.clone(), // each thread gets a clone
        move |pb, _| {
            let registry = Arc::clone(&registry);
            let mut timeline = Timeline::new(times, registry);

            let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
            let mut t_idx = 0;
            simulator.co_simulate(
                &mut rng(),
                sim_times,
                |t, tinc, _, w| {
                    while t_idx < times.len() && t + tinc >= times[t_idx] {
                        let structure = w.current_structure();
                        timeline.assign_structure(t_idx, &structure);
                        t_idx += 1;
                    }
                    true
                },
            );

            pb.inc(1);
            timeline
        },
    )

}



