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

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_eval_input(&cli.eval.input, is_rna)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }

    //NOTE: we use MPT as it is the more general method, but if you do a lot of
    //single-stranded evaluations, then it is a bit of overhead compared to
    //PairTable. May be worth to refactor and wrap PT and MPT into an enum which
    //implements LoopDecomposition.
    let pairings = MultiPairTable::try_from(&structure)?;
    let energy = model.energy_of_structure(&sequence, &pairings)?;

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

