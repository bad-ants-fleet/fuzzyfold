use std::io::Write;
use log::info;
use colored::*;
use env_logger::Builder;
use clap::Args;
use clap::Parser;
use clap::ArgAction;
use anyhow::Result;
use ahash::AHashMap;

use fuzzyfold::energy::EnergyModel;
use fuzzyfold::structure::DotBracketVec;
use fuzzyfold::structure::PairTable;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics::LoopStructure;
use fuzzyfold::kinetics::reaction::ApplyMove;


#[derive(Debug, Args)]
pub struct LocminInput {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// Input file (FASTA-like), or "-" for stdin
    #[arg(short, long, default_value_t = 0.0)]
    pub delta: f64,

    /// Verbosity (-v = info, -vv = debug)
    #[arg(short, long, action = ArgAction::Count)]
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

fn find_neighbors<'a, E: EnergyModel>(
    lss: &mut LoopStructure<'a, E>,
    structure: &DotBracketVec,
    max_delta: i32,
    neighbors: &mut AHashMap<DotBracketVec, i32>,
) {
    let mut mdbr = structure.clone();
    for (bp_move, delta) in lss.all_moves() {
        mdbr.apply_move(bp_move);
        if neighbors.contains_key(&mdbr) {
            mdbr.undo_move(bp_move);
            continue;
        } 

        lss.apply_move(bp_move);
        neighbors.insert(mdbr.clone(), lss.energy());

        if delta <= max_delta {
            println!("{} {:>6.2}", lss, (lss.energy() as f64 / 100.0));
            find_neighbors(lss, structure, max_delta - delta, neighbors);
        } else {
            //println!("{} {:>6.2} *", lss, (lss.energy() as f64 / 100.0));
        }
        lss.undo_move(bp_move);
        mdbr.undo_move(bp_move);
    }
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    init_logging(cli.lmin.verbose);

    let emodel = cli.energy.build_model();

    let (header, sequence, structure) = read_eval_input(&cli.lmin.input)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    } else {
        println!(">LM delta={:<.2}", cli.lmin.delta)
    }
    println!("{}", sequence);


    let pairings = PairTable::try_from(&structure)?;
    let mut lss = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
    println!("{} {:>6.2}", lss, (lss.energy() as f64 / 100.0));
    let mut neighbors = AHashMap::default();
    let delta = (cli.lmin.delta * 100.).round() as i32;
    find_neighbors(&mut lss, &structure, delta, &mut neighbors);

    Ok(())

}

