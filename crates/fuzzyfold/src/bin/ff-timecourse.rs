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
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_kinetics::Metropolis;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::timeline::Timeline;
use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistry;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelParams;
use fuzzyfold::kinetics_parsers::TimelineParameters;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
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
    kinetics: RateModelParams,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0);

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;

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

    let times = cli.simulation.get_output_times();
    let mut registry = MacrostateRegistry::from((&sequence, &emodel));
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

                let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
                let mut simulator = LoopStructureSSA::from((loops, &rmodel));
                let mut t_idx = 0;
                simulator.simulate(
                    &mut rng(),
                    cli.simulation.t_end,
                    |t, tinc, _, ls| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            let structure = DotBracketVec::from(ls);
                            timeline.assign_structure(t_idx, &structure);
                            t_idx += 1;
                        }
                        true
                    },
                );

                pb.inc(1);
                timeline
            },
        ).collect();
    pb.finish_with_message("All simulations complete!");

    // Master timeline
    for timeline in timelines {
        master.merge(timeline);
    }

    println!("Final Timeline:\n{}", master);
    plot_occupancy_over_time(&master, &format!("ff_{}.svg", name), cli.simulation.t_ext, cli.simulation.t_end);

    if let Some(path) = cli.timeline {
        let serial = master.to_serializable();
        let json = to_string_pretty(&serial)?;
        fs::write(path, json)?;
    }

    Ok(())
}



