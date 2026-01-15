use rand::rng;
use clap::Parser;
use anyhow::Result;
use colored::*;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::LoopStructure;
use ff_kinetics::LoopStructureSSA;

use fuzzyfold::input_parsers::read_fasta_like_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
//TODO: support seeded rng.

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

    /// Set the number of simulation steps (ignore t_end!!)
    #[arg(long, default_value_t = 0)]
    num_steps: usize,

    /// Do not print trajectory, only last structure.
    #[arg(short, long)]
    pub silent: bool, 

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = cli.energy.build_model();
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;
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

    let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
    let mut simulator = LoopStructureSSA::from((loops, &rmodel));

    if cli.silent { // Mostly for benchmarking anyway.
        if cli.num_steps == 0 { 
            simulator.simulate(&mut rng(), cli.t_end,
                |_, _, _, _|  true);
            let ls = simulator.current_structure();
            println!("{} {:8.2}", ls, ls.energy() as f64 / 100.);
        } else {
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, |_, _, _, _| {
                steps += 1;
                steps < cli.num_steps 
            });
            let ls = simulator.current_structure();
            println!("{} {:8.2}", ls, ls.energy() as f64 / 100.);
        }
    } else if cli.num_steps == 0 { 
        let callback = |t, tinc, flux, ls: &LoopStructure<'_, _>| {
            println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                ls,
                ls.energy() as f64 / 100.,
                t,
                tinc,
                1.0 / flux,
            );
            true
        };
        simulator.simulate(&mut rng(), cli.t_end, callback);
    } else {
        let mut steps = 0;
        let callback = |t, tinc, flux, ls: &LoopStructure<'_, _>| {
            println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                ls,
                ls.energy() as f64 / 100.,
                t,
                tinc,
                1.0 / flux,
            );
            steps += 1;
            steps < cli.num_steps 
        };
        simulator.simulate(&mut rng(), f64::MAX, callback);
    }

    Ok(())
}



