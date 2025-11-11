use clap::Args;
use clap::Parser;
use colored::*;
use anyhow::Result;

use rand::rng;
use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;
use ff_kinetics::Metropolis;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
//TODO: support seeded rng.

#[derive(Debug, Args)]
pub struct RateModelParams {
    /// Metropolis rate constant (must be > 0).
    #[arg(long, default_value_t = 1.0)]
    pub k0: f64,
}

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelParams,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0);

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }
    println!("{} {:>8} {:>14} {:>14} {:>15}",
        sequence,
        "energy".green(),
        "arivaltime".cyan(),
        "waitingtime".cyan(),
        "mean-waiting".cyan(),
    );

    let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, &rmodel));

    simulator.simulate(
        &mut rng(), 
        cli.t_end, 
        |t, tinc, flux, ls| {
            println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                ls,
                ls.energy() as f64 / 100.,
                t,
                tinc,
                1.0 / flux,
            );
            true
        },
    );
    Ok(())
}



