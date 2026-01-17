use ff_structure::DotBracketVec;
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
#[command(version, about = "Stochastic Simulation Algorithm for cotranscriptional RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Extension time
    #[arg(long, default_value_t = 0.2)]
    t_ext: f64,

    ///Pausing sites 
    #[arg(long, value_delimiter = ',')]
    t_p: Option<Vec<f64>>,

    ///Pausing sites positions  
    #[arg(long, value_delimiter = ',')]
    p_pos: Option<Vec<usize>>,

    ///Postranscriptional folding time 
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

    /// Set the number of simulation steps (ignore t_end!!)
    #[arg(long, default_value_t = 0)]
    num_steps: usize,

    /// Input structure 
    #[arg(long, default_value = ".")]
    start_structure: String,


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
    let pairings = PairTable::try_from(cli.start_structure.as_str())?;
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


    //build times vector 
    let mut times: Vec<f64> = Vec::new();
    let start = cli.start_structure.len();
    let mut idx = 0;

    if let Some(pause_times) = &cli.t_p {
       if let Some(pause_positions) = &cli.p_pos {
            let mut pos = start;
            times.push(cli.t_ext);
            pos += 1;
            for &p in pause_positions {
                while pos < p {
                    times.push(times.last().unwrap() + cli.t_ext);
                    pos += 1;
                }
                if pos == p {
                    times.push(times.last().unwrap() + pause_times[idx]);
                    idx += 1;
                    pos += 1;
                }
            }
            while pos < (sequence.len() -1) {
                times.push(times.last().unwrap() + cli.t_ext);
                pos += 1;
            }
            if pos == (sequence.len() - 1) {
                times.push(times.last().unwrap() + cli.t_end);
            }
       }
    } else {

        times.push(cli.t_ext);
        for _ in (start + 1)..(sequence.len() - 1) {
            times.push(times.last().unwrap() + cli.t_ext);
        }
        times.push(times.last().unwrap() + cli.t_end);
    }

    let mut steps = 0;

    simulator.co_simulate(
        &mut rng(), 
        times, 
        |t, tinc, flux, ls| {
            println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                ls,
                ls.energy() as f64 / 100.,
                t,
                tinc,
                1.0 / flux,
            );
            steps += 1;
            if steps == cli.num_steps {
               return false;
            }
            true
        },
    );
    Ok(())
}

