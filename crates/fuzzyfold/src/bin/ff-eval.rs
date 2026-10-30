use std::io::Write;
use log::{info, log_enabled, Level};
use colored::*;
use clap::Args;
use clap::Parser;
use clap::ArgAction;
use anyhow::Result;
use env_logger::Builder;

use fuzzyfold::energy::EnergyModel;
use fuzzyfold::structure::MultiPairTable;
use fuzzyfold::structure::ClosedRegions;
use fuzzyfold::input_parsers::ruler;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::{log_terms, EnergyModelArguments};


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
#[command(name = "ff-eval", version, about = "The free energy of a structure, also with pseudoknots and several strands.", after_help = fuzzyfold::CITATION)]
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
    fuzzyfold::restore_sigpipe();
    let cli = Cli::parse();
    init_logging(cli.eval.verbose);

    let model = cli.energy.build_model()?;

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_eval_input(&cli.eval.input, is_rna)?;
    cli.energy.warn_circular(&sequence);
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    let pairings = MultiPairTable::try_from(&structure)?;
    if log_enabled!(Level::Info) {
        log_terms(&model.energy_terms(&sequence, &pairings)?);
    }
    let energy = model.energy_of_structure(&sequence, &pairings)?;

    if log_enabled!(Level::Info) && !pairings.is_nested() {
        info!("Pseudoknot density: {}", ClosedRegions::from(&pairings).density());
    }
    info!("{}", ruler(sequence.len() - 1).magenta());
    println!("{}\n{} {}", sequence, structure, format!("{:>6.2}", energy as f64 / 100.0).green());
    info!("{}", ruler(sequence.len() - 1).magenta());

    if sequence.has_indistinguishable_strands() {
        // St1&St2 vs St2&St1
        // AAA&AAA vs AAA&AAA
        // ...&... vs ...&... -> 1x
        // .(.&.). vs .(.&.). -> 1x
        // (.)&... vs ...&(.) -> 2x
        info!("{}", "Note: Due to physically indistinguishable input sequences, the 
ensemble frequencies of distinguishable conformations do not 
correspond to the frequencies of dot-bracket representations. 
Symmetric structures have only one representation, while 
asymmetic structures have R representations.".red().bold())
    }
    Ok(())

}

