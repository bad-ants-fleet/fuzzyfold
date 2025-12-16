use std::io::Write;
use log::info;
use colored::*;
use env_logger::Builder;
use clap::Args;
use clap::Parser;
use clap::ArgAction;
use anyhow::Result;

use fuzzyfold::energy::EnergyModel;
use fuzzyfold::structure::MultiPairTable;
use fuzzyfold::input_parsers::ruler;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;


#[derive(Debug, Args)]
pub struct EvalInput {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    pub input: String,

    /// Verbosity (-v = info, -vv = debug)
    #[arg(short, long, action = ArgAction::Count)]
    pub verbose: u8,
}


#[derive(Debug, Parser)]
#[command(name = "ff-eval")]
#[command(author, version, about)]
pub struct Cli {
    #[command(flatten)]
    pub eval: EvalInput,

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
    init_logging(cli.eval.verbose);

    let model = cli.energy.build_model();

    let (header, sequence, structure) = read_eval_input(&cli.eval.input)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let pairings = MultiPairTable::try_from(&structure)?;
    let energy = model.energy_of_structure(&sequence, &pairings);

    info!("{}", ruler(sequence.len() - 1).magenta());
    println!("{}\n{} {}", sequence, structure, format!("{:>6.2}", energy as f64 / 100.0).green());
    info!("{}", ruler(sequence.len() - 1).magenta());

    Ok(())
}

