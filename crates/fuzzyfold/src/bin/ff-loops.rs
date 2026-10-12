use std::io::{self, BufRead};

use anyhow::{Result, anyhow};
use clap::Parser;

use fuzzyfold::energy::parse_structure;
use fuzzyfold::energy::{
    density_annotation, extended_dot_bracket_to_pair_table, is_complex, rastegari_density,
    zhao_density,
};
use ff_structure::ExtendedDotBracketVec;

#[derive(Debug, Parser)]
#[command(name = "ff-loops")]
#[command(about = "Decompose an extended dot-bracket structure into its loop types.\n\
                   Prints one loop per line: type, location status, closing pairs,\n\
                   and pseudoloop-specific fields (bands, pup, nested).")]
struct Cli {
    /// Extended dot-bracket string, or \"-\" to read from stdin
    #[arg(value_name = "STRUCTURE", default_value = "-")]
    structure: String,

    /// Report the density of each pseudoloop and of the structure
    #[arg(short, long)]
    density: bool,
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    let cli = Cli::parse();

    let st_str = if cli.structure == "-" {
        let stdin = io::stdin();
        let mut lines = stdin.lock().lines();
        loop {
            match lines.next() {
                None => return Err(anyhow!("No structure on stdin")),
                Some(l) => {
                    let l = l?.trim().to_string();
                    if l.is_empty() || l.starts_with('>') {
                        continue;
                    }
                    // skip a bare sequence line (only letters, no brackets)
                    if l.chars().all(|c| c.is_ascii_alphabetic()) {
                        continue;
                    }
                    break l;
                }
            }
        }
    } else {
        cli.structure.clone()
    };

    let loops = parse_structure(&st_str)?;

    for (i, lp) in loops.iter().enumerate() {
        println!("{i:3}  {lp}");
    }

    if cli.density {
        let edb = ExtendedDotBracketVec::try_from(st_str.as_str())?;
        let pt = extended_dot_bracket_to_pair_table(&edb)?;

        println!("{}", "-".repeat(60));
        for d in density_annotation(&pt) {
            println!(
                "region {:>4}-{:<4}  bands {}  zhao {}  rastegari {}  band_lengths {:?}",
                d.region.0, d.region.1, d.n_bands, d.zhao, d.rastegari, d.band_lengths
            );
        }
        println!(
            "density  zhao {}  rastegari {}{}",
            zhao_density(&pt),
            rastegari_density(&pt),
            if is_complex(&pt) {
                "  complex (outside the density-2 class)"
            } else {
                ""
            }
        );
    }

    Ok(())
}
