use std::io::{stdin, Write, BufRead};
use env_logger::Builder;
use clap::{Parser, Args, ArgAction};
use anyhow::Result;
use anyhow::anyhow;

use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_structure::PairTable;
use ff_structure::DotBracketVec;
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

    let mut lines = stdin().lock().lines();

    // --- First line: sequence and maybe extra numbers ---
    let first = lines.next().ok_or_else(|| anyhow!("Missing first line"))??;

    let mut parts = first.split_whitespace();
    let token = parts.next().ok_or_else(|| anyhow!("Missing sequence"))?;
    let sequence = NucleotideVec::from_lossy(token);
    let _: Vec<f64> = parts.map(|s| s.parse().unwrap()).collect();
    println!("{}", sequence);

    // --- Remaining lines: dot-bracket + energy ---
    for line in lines {
        let line = line?;
        let mut parts = line.split_whitespace();

        let token = parts.next().ok_or_else(|| anyhow!("Missing structure"))?;
        let structure = DotBracketVec::try_from(token)?;
        let ref_en: f64 = parts.next().ok_or_else(|| anyhow!("Missing energy"))?.parse()?;

        let pairings = PairTable::try_from(&structure)?;
        let energy = model.energy_of_structure(&sequence, &pairings);

        let mark = if (ref_en * 100f64).round() as i32 != energy { "*" } else { "" };
        if mark == "*" {
            println!("{} {:6.2} {:6.2} {}", structure, energy as f64 / 100.0, ref_en, mark);
        }
    }

    Ok(())
}


