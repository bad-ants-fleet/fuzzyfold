use std::io::Write;
use colored::*;
use env_logger::Builder;
use clap::Args;
use clap::Parser;
use clap::ArgAction;
use clap::ArgGroup;
use anyhow::Result;

use fuzzyfold::structure::PairTable;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics::LoopStructure;


#[derive(Debug, Args)]
#[command(
    group = ArgGroup::new("criterion")
        .required(true)
        .multiple(true)
        .args(["delta", "maxdist"])
)]
pub struct LocminInput {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// give delta energy
    #[arg(short, long)]
    pub delta: Option<f64>,

    /// Input file (FASTA-like), or "-" for stdin
    #[arg(short, long)]
    pub maxdist: Option<usize>,

    /// give delta energy
    #[arg(short, long)]
    pub sorted: bool, 

    /// Verbosity (-v = info, -vv = debug)
    #[arg(short, long, hide = true, action = ArgAction::Count)]
    pub verbose: u8,
}


#[derive(Debug, Parser)]
#[command(name = "ff-eval")]
#[command(author, version, about)]
pub struct Cli {
    #[command(flatten)]
    pub lmin: LocminInput,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    pub energy: EnergyModelArguments,
}

fn init_logging(verbosity: u8) {
    let level = match verbosity {
        0 => "warn",
        1 => "info",
        _ => "debug",
    };

    Builder::from_env(env_logger::Env::default().default_filter_or(level))
        .format(|buf, record| {
            // no prefix, just the message
            writeln!(buf, "{}", record.args())
        })
        .init();
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.lmin.verbose);
    let emodel = cli.energy.build_model();
    let (header, sequence, structure) = read_eval_input(&cli.lmin.input)?;
    let pairings = PairTable::try_from(&structure)?;

    let (delta, distance, info) = match (cli.lmin.delta, cli.lmin.maxdist) {
        (Some(d), None) => {
            ((d * 100.0) as i32, usize::MAX, format!("delta = {:<.2}", d))
        },
        (None, Some(n)) => (i32::MAX/2, n, format!("maxdist = {}", n)),
        (Some(d), Some(n)) => ((d * 100.0) as i32, n, format!("delta = {:<.2} maxdist = {}", d, n)),
        _ => unreachable!("clap guarantees one is set"),
    };

    if let Some(h) = header {
        println!("{} ({})", h.yellow(), info);
        println!("{}", sequence);
    } else {
        println!(">LM ({})", info);
        println!("{}", sequence);
    }

    let mut lss = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();

    if !cli.lmin.sorted {
        lss.generate_neighbors(delta, distance, 
            |db, en| { 
                println!("{} {:.2}", db, en as f64 / 100.0);
            }
        );
    } else {
        let mut neighbors = Vec::new();
        lss.generate_neighbors(delta, distance, 
            |db, en| { 
                neighbors.push((db.clone(), en)); }
        );
        neighbors.sort_by_key(|(_, en)| *en);
        for (db, en) in &neighbors {
            println!("{} {:.2}", db, *en as f64 / 100.0);
        }
    }

    Ok(())

}

