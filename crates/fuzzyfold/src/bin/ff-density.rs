//! Report the density of an extended dot-bracket structure.
//!
//! Density is Zhao's (2005) eq. 1: the largest number of bands of a pseudoloop
//! crossing one position, maximised over pseudoloops. It is not the band count
//! -- a kissing hairpin has three bands and density 2.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};

use anyhow::{Context, Result};
use clap::Parser;

use ff_structure::ExtendedDotBracketVec;
use fuzzyfold::energy::{
    density_annotation, extended_dot_bracket_to_pair_table, is_complex, rastegari_density,
    zhao_density,
};

#[derive(Debug, Parser)]
#[command(name = "ff-density")]
#[command(about = "Report the density of an extended dot-bracket structure.\n\
                   Prints one Zhao density per input structure; --tsv adds the \
                   other measures. Reads one structure per line, so it can \
                   annotate a whole file.")]
struct Cli {
    /// Extended dot-bracket string, or "-" to read one structure per line from stdin
    #[arg(value_name = "STRUCTURE", default_value = "-")]
    structure: String,

    /// Read structures from this file instead, one per line
    #[arg(short, long, value_name = "PATH")]
    file: Option<String>,

    /// Print a header and all measures instead of the density alone
    #[arg(short, long)]
    tsv: bool,

    /// One row per pseudoloop rather than one per structure
    #[arg(short, long)]
    per_pseudoloop: bool,
}

/// Columns are stable: downstream tables join on them.
const STRUCTURE_HEADER: &str = "len\tpairs\tpseudoloops\tbands\tzhao\trastegari\tcomplex";
const PSEUDOLOOP_HEADER: &str = "struct\ti\tj\tbands\tzhao\trastegari\tband_lengths";

/// Skips what a .dp or FASTA-like file puts around the structures: blank lines,
/// `#` comments, `>` headers, and sequence lines. Letters are bracket families
/// now, so a sequence is only recognisable by having no other bracket.
fn is_structure(line: &str) -> bool {
    !line.is_empty()
        && !line.starts_with('#')
        && !line.starts_with('>')
        && line.chars().any(|c| ".()[]{}<>".contains(c))
}

/// `Ok(None)` on a structure this crate cannot parse: batch runs report it and
/// carry on rather than abandoning the remaining lines.
fn report(out: &mut impl Write, cli: &Cli, nth: usize, s: &str) -> Result<Option<()>> {
    let pt = match ExtendedDotBracketVec::try_from(s).and_then(|edb| {
        extended_dot_bracket_to_pair_table(&edb)
    }) {
        Ok(pt) => pt,
        Err(e) => {
            eprintln!("structure {nth}: {e}");
            writeln!(out, "NA")?;
            return Ok(None);
        }
    };

    let annotation = density_annotation(&pt);

    if cli.per_pseudoloop {
        for d in &annotation {
            let lengths: Vec<String> = d.band_lengths.iter().map(|l| l.to_string()).collect();
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                nth,
                d.region.0,
                d.region.1,
                d.n_bands,
                d.zhao,
                d.rastegari,
                lengths.join(",")
            )?;
        }
    } else if cli.tsv {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pt.len(),
            pt.iter().filter(|p| p.is_some()).count() / 2,
            annotation.len(),
            annotation.iter().map(|d| d.n_bands).sum::<usize>(),
            zhao_density(&pt),
            rastegari_density(&pt),
            if is_complex(&pt) { "yes" } else { "no" }
        )?;
    } else {
        writeln!(out, "{}", zhao_density(&pt))?;
    }

    Ok(Some(()))
}

/// One structure per line, skipping whatever else the file carries.
fn report_stream(out: &mut impl Write, cli: &Cli, reader: impl BufRead) -> Result<bool> {
    let mut failed = false;
    let mut nth = 0;
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if !is_structure(line) {
            continue;
        }
        nth += 1;
        if report(out, cli, nth, line)?.is_none() {
            failed = true;
        }
    }
    Ok(failed)
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    if cli.tsv {
        writeln!(out, "{STRUCTURE_HEADER}")?;
    } else if cli.per_pseudoloop {
        writeln!(out, "{PSEUDOLOOP_HEADER}")?;
    }

    let failed = match cli.file.as_deref() {
        Some(path) => {
            let f = File::open(path).with_context(|| format!("cannot open '{path}'"))?;
            report_stream(&mut out, &cli, BufReader::new(f))?
        }
        None if cli.structure == "-" => report_stream(&mut out, &cli, io::stdin().lock())?,
        None => report(&mut out, &cli, 1, &cli.structure)?.is_none(),
    };

    out.flush()?;
    if failed {
        std::process::exit(1);
    }
    Ok(())
}
