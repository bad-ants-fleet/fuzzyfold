use std::fs::File;
use std::io::{BufRead, BufReader, IsTerminal, Write, stdin, stdout};
use std::path::Path;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use clap::{Parser, ValueEnum};
use resvg::usvg::fontdb;
use resvg::{tiny_skia, usvg};

use ff_structure::PairTable;
use ff_structure::arc_diagram::arc_diagram_svg_opt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Format {
    /// Rasterized bitmap (default)
    Png,
    /// Vector output
    Svg,
}

#[derive(Debug, Parser)]
#[command(name = "ff-arc-diagram")]
#[command(
    about = "Render a secondary structure as an arc diagram (PNG or SVG).",
    long_about = "Render a secondary structure as an arc diagram.\n\
\n\
Input is FASTA-like: an optional >header, an optional sequence line, and a\n\
dot-bracket structure line. A lone structure line is accepted — bases are then\n\
drawn as neutral dots without nucleotide labels.\n\
\n\
Extended dot-bracket is supported, so pseudoknots render as crossing arcs:\n\
([{<ABCD and their closing partners.\n\
\n\
Output is PNG by default. The format follows the --output extension, and\n\
--format overrides it (needed when writing to stdout).\n\
\n\
Examples:\n\
  ff-arc-diagram structure.fa                      # -> arc.png\n\
  ff-arc-diagram structure.fa -o fig.svg           # vector, from the extension\n\
  echo '((([[[)))...]]]' | ff-arc-diagram -o pk.png\n\
  echo '((([[[)))...]]]' | ff-arc-diagram --scale 4\n\
  ff-arc-diagram in.fa --format svg -o - > arc.svg"
)]
struct Cli {
    /// Input file, or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Output path, or "-" for stdout
    #[arg(short, long, value_name = "OUTPUT", default_value = "arc.png")]
    output: String,

    /// Output format; defaults to the --output extension, else png
    #[arg(short, long, value_enum)]
    format: Option<Format>,

    /// Pixels per SVG unit when rasterizing (png only)
    #[arg(long, default_value = "2.0", value_name = "FACTOR")]
    scale: f32,
}

/// A line is a structure if it holds only dot-bracket characters and at least
/// one of them is not a letter; anything else is taken as the sequence. A line
/// of letters alone is indistinguishable from a sequence, and family 0 is
/// always `()`, so a real structure always carries a non-letter bracket.
fn is_structure(line: &str) -> bool {
    !line.is_empty()
        && line
            .chars()
            .all(|c| ".()[]{}<>+&".contains(c) || c.is_ascii_alphabetic())
        && line.chars().any(|c| ".()[]{}<>".contains(c))
}

/// Reads an optional sequence and a required structure. Accepts either order,
/// a lone structure, and skips >headers and blank lines.
fn read_input(path: &str) -> Result<(Option<String>, String)> {
    let reader: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        Box::new(BufReader::new(
            File::open(path).with_context(|| format!("cannot open '{path}'"))?,
        ))
    };

    let mut sequence = None;
    let mut structure = None;

    for line in reader.lines() {
        let line = line?.trim().replace(' ', "");
        if line.is_empty() || line.starts_with('>') {
            continue;
        }
        if is_structure(&line) {
            if structure.is_none() {
                structure = Some(line);
            }
        } else if sequence.is_none() {
            sequence = Some(line);
        }
        if structure.is_some() && sequence.is_some() {
            break;
        }
    }

    let structure = structure.ok_or_else(|| {
        let src = if path == "-" { "stdin" } else { path };
        anyhow!("no dot-bracket structure line found in {src}")
    })?;
    Ok((sequence, structure))
}

/// Resolves the generic `monospace` family to a face that actually exists, so
/// the nucleotide labels do not silently fall back to a default sans face.
fn pick_monospace(db: &fontdb::Database) -> Option<String> {
    const CANDIDATES: &[&str] = &[
        "Menlo",
        "DejaVu Sans Mono",
        "Liberation Mono",
        "Noto Sans Mono",
        "Consolas",
        "Courier New",
    ];
    let has = |want: &str| {
        db.faces()
            .any(|f| f.families.iter().any(|(name, _)| name == want))
    };
    CANDIDATES
        .iter()
        .find(|c| has(c))
        .map(|c| c.to_string())
        .or_else(|| {
            db.faces()
                .find(|f| f.monospaced)
                .and_then(|f| f.families.first().map(|(name, _)| name.clone()))
        })
}

fn rasterize(svg: &str, scale: f32) -> Result<Vec<u8>> {
    if !(scale.is_finite() && scale > 0.0) {
        return Err(anyhow!("--scale must be a positive number, got {scale}"));
    }

    let mut db = fontdb::Database::new();
    db.load_system_fonts();
    if let Some(mono) = pick_monospace(&db) {
        db.set_monospace_family(mono);
    }

    let options = usvg::Options {
        fontdb: Arc::new(db),
        ..Default::default()
    };
    let tree = usvg::Tree::from_str(svg, &options).context("rendering the arc diagram failed")?;

    let size = tree.size();
    let w = (size.width() * scale).ceil() as u32;
    let h = (size.height() * scale).ceil() as u32;
    let mut pixmap = tiny_skia::Pixmap::new(w.max(1), h.max(1))
        .ok_or_else(|| anyhow!("{w}x{h} px is too large to allocate; lower --scale"))?;

    resvg::render(
        &tree,
        tiny_skia::Transform::from_scale(scale, scale),
        &mut pixmap.as_mut(),
    );
    Ok(pixmap.encode_png()?)
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let format = cli.format.unwrap_or_else(|| {
        match Path::new(&cli.output)
            .extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("svg") => Format::Svg,
            _ => Format::Png,
        }
    });

    let (sequence, structure) = read_input(&cli.input)?;

    if let Some(seq) = &sequence
        && seq.len() != structure.len()
    {
        return Err(anyhow!(
            "sequence is {} nt but structure is {} nt",
            seq.len(),
            structure.len()
        ));
    }

    let pt = PairTable::try_from(structure.as_str())?;
    let svg = arc_diagram_svg_opt(sequence.as_deref(), &pt);

    let bytes = match format {
        Format::Svg => svg.into_bytes(),
        Format::Png => rasterize(&svg, cli.scale)?,
    };

    if cli.output == "-" {
        if format == Format::Png && stdout().is_terminal() {
            return Err(anyhow!(
                "refusing to write PNG bytes to a terminal; redirect stdout or pass --format svg"
            ));
        }
        stdout().write_all(&bytes)?;
    } else {
        std::fs::write(&cli.output, &bytes)
            .with_context(|| format!("cannot write '{}'", cli.output))?;
        eprintln!(
            "Arc diagram written to '{}' ({:?}, {} bytes)",
            cli.output,
            format,
            bytes.len()
        );
    }
    Ok(())
}
