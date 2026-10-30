use clap::Parser;
use anyhow::Result;
use colored::*;
use rand::rng;
use std::sync::Arc;

use ff_structure::MultiPairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Elongation;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::SSA;
use ff_kinetics::extension_schedule;

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::input_parsers::read_fasta_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
//TODO: support seeded rng.

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Read input as standard FASTA (DNA alphabet, multi-line sequence, no structure)
    #[arg(long)]
    fasta: bool,

    /// Simulation time at the full-length molecule.
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

    /// Optional simulation time per nucleotide during co-transcriptional folding.
    #[arg(long)]
    t_ext: Option<f64>,

    /// Set the number of simulation steps (ignore t_end!!)
    #[arg(long, hide = true, default_value_t = 0)]
    num_steps: usize,

    /// Do not print trajectory, only last structure.
    #[arg(long)]
    silent: bool, 

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model()?);
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) =
        if cli.fasta {
            read_fasta_input(&cli.input, true)?  // is_rna is ignored, assuming DNA alphabet
        } else if cli.t_ext.is_some() {
            read_cotr_input(&cli.input, is_rna)?
        } else {
            match read_eval_input(&cli.input, is_rna) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };
    cli.energy.warn_circular(&sequence);

    let pairings = MultiPairTable::try_from(&structure)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }
    println!("{} {:>8} {:>14} {:>14} {:>15}",
        sequence,
        "energy".green(),
        "arrival-time".cyan(),
        "waiting-time".cyan(),
        "mean-waiting".cyan(),
    );

    let num_ext = sequence.len() - structure.len();
    let end = cli.t_ext.map_or(0.0, |t| t * num_ext as f64) + cli.t_end;
    let times = extension_schedule(cli.t_ext, num_ext, end).map_err(anyhow::Error::msg)?;

    let moves = LoopNeighbors::try_from((sequence, &pairings, emodel, rmodel.shift_moves()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    run_simulator(moves, rmodel, &times, cli.silent, cli.num_steps);

    Ok(())
}


fn run_simulator<W: Elongation + std::fmt::Display, K: RateModel>(
    moves: W,
    rmodel: K,
    times: &[f64],
    silent: bool,
    num_steps: usize,
)
{
    let mut simulator = SSA::from((moves, rmodel));

    match (silent, num_steps) {
        (true, 0) => {
            simulator.co_simulate(&mut rng(), times, |_, _, _, _|  true);
            let cs = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", cs, en as f64 / 100.);
        },
        (false, 0) => {
            simulator.co_simulate(&mut rng(), times, 
                |t, tinc, flux, w| {
                    println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                        w,
                        w.current_energy() as f64 / 100.,
                        t,
                        tinc,
                        1.0 / flux,
                    );
                    true
                });
        },
        (true, _) => {
            if times.len() > 1 {
                panic!("no support for co-transcriptional num-steps.")
            }
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, |_, _, _, _| {
                steps += 1;
                steps < num_steps 
            });
            let ls = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", ls, en as f64 / 100.);
        },
        (false, _) => {
            if times.len() > 1 {
                panic!("no support for co-transcriptional num-steps.")
            }
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, 
                |t, tinc, flux, w| {
                    println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                        w,
                        w.current_energy() as f64 / 100.,
                        t,
                        tinc,
                        1.0 / flux,
                    );
                    steps += 1;
                    steps < num_steps 
                });
        },
    }
}
