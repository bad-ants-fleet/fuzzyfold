//! `ff-pseudo-findpath` — PK-aware beam-search path between two RNA structures.
//!
//! Reads sequence + optional start + target from a file or stdin.
//! File format (whitespace-trimmed, empty lines ignored):
//!
//!   Line 1: RNA sequence (ACGU)
//!   Line 2: start structure in extended dot-bracket (or target if only 2 lines)
//!   Line 3: target structure in extended dot-bracket
//!
//! If only two non-empty lines are provided the start defaults to the fully-unpaired structure.

use ff_findpath::pseudo_findpath::findpath_pseudo;
use ff_energy::{ViennaRNA, parameters::{RNA_MT09, RNA_DP09}};

use std::io::{self, BufRead};
use std::fs;
use std::path::Path;

use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "ff-pseudo-findpath",
    about = "PK-aware beam search from a start to a pseudoknotted target structure",
    long_about = None,
)]
struct Args {
    /// Input file (sequence + optional start + target, one per line).
    /// If omitted, reads from stdin.
    #[arg(short, long, value_name = "FILE")]
    file: Option<String>,

    /// Search (beam) width. 1 = greedy; larger values explore more paths.
    #[arg(short = 'm', long = "beam-width", value_name = "N", default_value = "10")]
    beam_width: usize,

    /// Energy ceiling in kcal/mol. Intermediates above this are pruned.
    /// Defaults to no ceiling (unlimited).
    #[arg(long, value_name = "KCAL")]
    max_energy: Option<f64>,

    /// Temperature in Celsius for nearest-neighbor stacking parameters.
    #[arg(long, value_name = "CELSIUS", default_value = "37")]
    celsius: f64,

    /// Print per-step cache hit statistics to stderr.
    #[arg(long)]
    cache_stats: bool,
}

fn read_lines(source: &str) -> io::Result<Vec<String>> {
    if source == "-" {
        let stdin = io::stdin();
        Ok(stdin.lock().lines()
            .map(|l| l.map(|s| s.trim().to_owned()))
            .filter(|l| l.as_ref().map_or(true, |s| !s.is_empty()))
            .collect::<io::Result<Vec<_>>>()?)
    } else {
        let content = fs::read_to_string(Path::new(source))?;
        Ok(content.lines()
            .map(|l| l.trim().to_owned())
            .filter(|l| !l.is_empty())
            .collect())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // ── read input ────────────────────────────────────────────────────────────
    let source = args.file.as_deref().unwrap_or("-");
    let lines = read_lines(source)?;

    let (sequence, start, target) = match lines.len() {
        2 => (lines[0].as_str(), None, lines[1].as_str()),
        3.. => (lines[0].as_str(), Some(lines[1].as_str()), lines[2].as_str()),
        _ => {
            eprintln!("Error: need at least 2 non-empty lines (sequence + target)");
            std::process::exit(1);
        }
    };

    // ── build model ───────────────────────────────────────────────────────────
    let model = ViennaRNA::from_andrunescu_params(&RNA_MT09)
        .with_pseudoknot_params(RNA_DP09);

    // ── run findpath ──────────────────────────────────────────────────────────
    let (path, stats) = findpath_pseudo(
        &model,
        sequence,
        start,
        target,
        args.beam_width,
        args.max_energy,
        args.cache_stats,
    )
    .unwrap_or_else(|e| {
        eprintln!("findpath_pseudo failed: {e}");
        std::process::exit(1);
    });

    // ── print results ─────────────────────────────────────────────────────────
    println!("# ff-pseudo-findpath  beam={}", args.beam_width);
    println!("# sequence  {sequence}");
    if let Some(s) = start {
        println!("# start     {s}");
    }
    println!("# target    {target}");
    println!("#");
    println!("# {:>3}  {:<width$}  {:>10}  move", "i", "structure", "kcal/mol", width = target.len());

    for step in &path {
        let mv_str = step.move_applied.as_ref()
            .map(|m| m.to_string())
            .unwrap_or_else(|| "start".to_owned());
        println!(
            "  {:>3}  {:<width$}  {:>10.2}  {}",
            step.step_index,
            step.structure,
            step.energy,
            mv_str,
            width = target.len(),
        );
    }

    println!("#");
    println!("# saddle    {:>8.2} kcal/mol", stats.saddle_energy);
    println!("# barrier   {:>8.2} kcal/mol", stats.barrier_energy);
    println!("# end       {:>8.2} kcal/mol", stats.end_energy);

    Ok(())
}
