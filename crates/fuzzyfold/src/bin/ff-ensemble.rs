use std::sync::Arc;
use std::path::PathBuf;

use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use fuzzyfold::structure::MultiPairTable;
use fuzzyfold::energy::EnergyModel;
use fuzzyfold::kinetics::Elongation;
use fuzzyfold::kinetics::LoopNeighbors;
use fuzzyfold::kinetics::SSA;
use fuzzyfold::kinetics::RateModel; 

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimeArguments;
use fuzzyfold::kinetics::OutputTimes;
use fuzzyfold::kinetics::ensemble::{Ensemble, State};

#[derive(Parser)]
#[command(version, about = "The structures of many trajectories at the output times (.drf).", after_help = fuzzyfold::CITATION)]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimeArguments,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    fuzzyfold::restore_sigpipe();
    let cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model()?);
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) =
        if cli.simulation.t_ext.is_some() {
            read_cotr_input(&cli.input, is_rna)?
        } else {
            match read_eval_input(&cli.input, is_rna) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };
    cli.energy.warn_circular(&sequence);
    let num_ext = sequence.len() - structure.len();
    let pairings = MultiPairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow());
    }
    println!("{}", sequence);
 
    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let output = cli.simulation.output_times(num_ext)?;
    let schedule = cli.simulation.schedule(num_ext, &output)?;
    let times = output.times();
    
    let edits = rmodel.move_types().a_to_i;
    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.move_types()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    let trajectories: Vec<Vec<State>> =
        run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64, &output).collect();

    let mut ensemble = Ensemble::new(times);
    for trajectory in trajectories {
        ensemble.add_trajectory(trajectory);
    }
    println!("{}", "Finished simulations!".red());

    // With A-to-I edits, the sequence of each structure is a column.
    let drf_path = cli.output.with_extension("drf");
    ensemble.write_drf(&drf_path, edits)?;
    println!("Wrote drf file: {}", format!("{}", drf_path.display()).green()); 

    Ok(())
}

fn run_timecourse<W, K>(
    moves: W,
    rmodel: K,
    schedule: &[f64],
    num_sims: u64,
    output: &OutputTimes,
) -> impl ParallelIterator<Item = Vec<State>>
where
    W: Elongation + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
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
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let trajectory = simulator.sample(&mut rng(), schedule, output,
                    |w| (w.current_sequence().clone(), w.current_structure(), w.current_energy()));
                pb.inc(1);
                trajectory
            },
        )
}
