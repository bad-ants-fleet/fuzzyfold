use std::io::Write;
use std::sync::Arc;
use colored::*;
use clap::Args;
use clap::Parser;
use clap::ArgAction;
use clap::ArgGroup;
use anyhow::Result;
use env_logger::Builder;

use ff_energy::EnergyModel;
use ff_kinetics::shift_policy::*;
use fuzzyfold::structure::PairTable;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics::LoopNeighbors;


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

    /// Specify energy delta for enumeration.
    #[arg(short, long)]
    pub delta: Option<f64>,

    /// Specify maximum distance to starting structure.
    #[arg(short, long)]
    pub maxdist: Option<usize>,

    /// Enable three-way shift moves.
    #[arg(long)]
    pub three_way_shifts: bool, 

    /// Enable four-way shift moves.
    #[arg(long)]
    pub four_way_shifts: bool, 

    /// Retrurn an energetically sorted list.
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
    let emodel = Arc::new(cli.energy.build_model());
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

    match (cli.lmin.three_way_shifts, cli.lmin.four_way_shifts) {
        (false, false) => {
            let mut moves = LoopNeighbors::try_from(
                (sequence, &pairings, emodel, NoShift)
            ).expect("failed to build loop table");
            do_enumerate(&mut moves, delta, distance, cli.lmin.sorted);
        }
        (true, false) => {
            let mut moves = LoopNeighbors::try_from(
                (sequence, &pairings, emodel, ThreeWayOnly)
            ).expect("failed to build loop table");
            do_enumerate(&mut moves, delta, distance, cli.lmin.sorted);
        }
        (false, true) => {
            let mut moves = LoopNeighbors::try_from(
                (sequence, &pairings, emodel, FourWayOnly)
            ).expect("failed to build loop table");
            do_enumerate(&mut moves, delta, distance, cli.lmin.sorted);
        }
        (true, true) => {
            let mut moves = LoopNeighbors::try_from(
                (sequence, &pairings, emodel, ThreeAndFour)
            ).expect("failed to build loop table");
            do_enumerate(&mut moves, delta, distance, cli.lmin.sorted);
        }
    }
    Ok(())
}

fn do_enumerate<E: EnergyModel, P: ShiftPolicy>(
    moves: &mut LoopNeighbors<E, P>,
    delta: i32,
    distance: usize,
    sorted: bool,
) {
    if !sorted {
        moves.generate_neighbors(delta, distance, 
            |db, en| { 
                println!("{} {:.2}", db, en as f64 / 100.0);
            }
        );
    } else {
        let mut neighbors = Vec::new();
        moves.generate_neighbors(delta, distance, 
            |db, en| { 
                neighbors.push((db.clone(), en)); }
        );
        neighbors.sort_by_key(|(_, en)| *en);
        for (db, en) in &neighbors {
            println!("{} {:.2}", db, *en as f64 / 100.0);
        }
    }
}


