use std::fs::File;
use std::io::{BufRead, BufReader, stdin};

use anyhow::{Result, anyhow};
use clap::{Parser, ValueEnum};

use fuzzyfold::energy::parameters::{
    RNA_ANDRONESCU_2007, RNA_DP03, RNA_DP09, RNA_EXTENDED, RNA_MT09, RNA_TURNER_2004,
};
use fuzzyfold::energy::{NucleotideVec, PseudoEnergyModel, ViennaRNA, parse_structure};

/// Preset (nearest-neighbour, pseudoknot) pairings.
///
/// Naming: lower case (`dp03`, `dp09`) is a set of the 11 Dirks-Pierce
/// parameters on its own, as passed to `--pk`. UPPER CASE is a combination of
/// a nearest-neighbour model with one of those sets, so a preset name reads
/// `<NN>_<pk>`, with `T04` = `turner2004` and `A07` = `andronescu2007` as
/// `--nn` spells them. `DP09` keeps the published name rather than spelling
/// out `MT09_dp09`.
///
/// The `--nn` value is named, not the model, because `turner2004` and
/// `extended` are both Turner 2004 — `extended` adds the abasic pair
/// extension — so the prose name alone does not identify a parameter set.
///
/// Presets are shorthands only: `--nn` and `--pk` override either half, so
/// every NN x PK combination stays reachable without one.
///
/// Only `DP09` reproduces a trained pairing (Andronescu et al. 2010, Table 3
/// row 9). The paper's other named set is DP03 = MT99 + dp03 (row 2), where
/// MT99 is Mathews et al. (1999). This crate carries no MT99 parameters, so
/// DP03 cannot be expressed here at all; every other preset is cross-model.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum PkChoice {
    /// --nn mt09 --pk dp09 — the trained pairing (default)
    #[value(name = "DP09")]
    Dp09,
    /// --nn turner2004 --pk dp03 — cross-model; published DP03 is MT99 + dp03
    #[value(name = "T04_dp03")]
    T04Dp03,
    /// --nn turner2004 --pk dp09 — cross-model
    #[value(name = "T04_dp09")]
    T04Dp09,
    /// --nn andronescu2007 --pk dp03 — cross-model
    #[value(name = "A07_dp03")]
    A07Dp03,
}

/// Nearest-neighbour parameter set.
///
/// DNA (Mathews 2004) is deliberately absent: the Dirks-Pierce pseudoknot
/// parameters were trained on RNA only, so no DNA + DP combination is
/// scientifically meaningful.
#[derive(Debug, Clone, Copy, ValueEnum)]
enum NnChoice {
    /// Andronescu 2009 (MT09) — trained with dp09 as DP09, fixed at 37 °C
    Mt09,
    /// Turner 2004 — temperature dependent via --celsius
    Turner2004,
    /// Andronescu 2007 — fixed at 37 °C
    Andronescu2007,
    /// Turner 2004 + abasic pair extension — temperature dependent
    Extended,
}

/// Dirks-Pierce pseudoknot parameter set (Andronescu et al. 2010, Table 6).
#[derive(Debug, Clone, Copy, ValueEnum)]
enum DpChoice {
    /// Retrained 2009 values — best performing (F-measure 79%)
    Dp09,
    /// Original 2003 values — Dirks & Pierce, as used by HotKnots / NUPACK
    Dp03,
}

#[derive(Debug, Parser)]
#[command(name = "ff-calc-pseudo")]
#[command(
    about = "Evaluate pseudoknot free energy using the Dirks-Pierce model.",
    long_about = "Evaluate pseudoknot free energy using the Dirks-Pierce model.\n\
\n\
Parameters are chosen either by preset or component-wise. Each preset is\n\
exactly the --nn/--pk pair shown beside it:\n\
\n\
DP09      =  --nn mt09            --pk dp09     (trained pairing, default)\n\
T04_dp03  =  --nn turner2004      --pk dp03\n\
T04_dp09  =  --nn turner2004      --pk dp09\n\
A07_dp03  =  --nn andronescu2007  --pk dp03\n\
\n\
--nn  mt09 | turner2004 | andronescu2007 | extended\n\
--pk  dp09 | dp03\n\
\n\
Lower case (dp03, dp09) is a set of the 11 Dirks-Pierce parameters on its own,\n\
which is what --pk takes. UPPER CASE is a combination of a nearest-neighbour\n\
model with one of them, so a preset reads <NN>_<pk>, with T04 = turner2004 and\n\
A07 = andronescu2007. DP09 keeps its published name instead of MT09_dp09.\n\
Note that turner2004 and extended are different models: extended is Turner\n\
2004 plus the abasic pair extension, so \"Turner 2004\" alone is ambiguous.\n\
\n\
--nn and --pk override the corresponding half of the preset, so all eight\n\
RNA combinations are reachable.\n\
\n\
Only MT09 + dp09 is a trained pairing: DP09, Andronescu et al. (2010) Table 3\n\
row 9. The paper's DP03 is MT99 + dp03 (row 2); MT99 (Mathews et al. 1999) is\n\
not carried here, so DP03 cannot be expressed by this crate. Turner 2004 is a\n\
different nearest-neighbour model, not a substitute for it -- comparing\n\
turner2004 + dp03 against a published DP03 number differs in the NN model as\n\
well. Every combination other than MT09 + dp09 mixes independently fitted sets."
)]
struct Cli {
    /// Input file (FASTA-like: optional >header, sequence, structure), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Preset (NN, PK) pairing
    #[arg(long, value_enum, default_value = "DP09")]
    pk_params: PkChoice,

    /// Override the nearest-neighbour model of the preset
    #[arg(long, value_enum)]
    nn: Option<NnChoice>,

    /// Override the pseudoknot parameter set of the preset
    #[arg(long, value_enum)]
    pk: Option<DpChoice>,

    /// Temperature in Celsius (only affects turner2004 / extended; MT09 and A07 are fixed at 37 °C)
    #[arg(long, default_value = "37.0")]
    celsius: f64,

    /// Print the resolved parameter combination to stderr
    #[arg(long)]
    show_params: bool,

    /// Print the energy contribution of each loop individually
    #[arg(short = 'v', long)]
    verbose: bool,
}

fn read_input(path: &str) -> Result<(String, String)> {
    let reader: Box<dyn BufRead> = if path == "-" {
        Box::new(BufReader::new(stdin()))
    } else {
        Box::new(BufReader::new(File::open(path)?))
    };

    let mut sequence: Option<String> = None;
    let mut structure: Option<String> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            continue;
        }
        if sequence.is_none() {
            sequence = Some(line.replace(' ', ""));
        } else {
            structure = Some(line.replace(' ', ""));
            break;
        }
    }

    let seq = sequence.ok_or_else(|| anyhow!("Missing sequence line"))?;
    let st = structure.ok_or_else(|| anyhow!("Missing structure line"))?;
    Ok((seq, st))
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    let cli = Cli::parse();

    let (seq_str, st_str) = read_input(&cli.input)?;

    let seq = NucleotideVec::try_from_rna(&seq_str)?;
    let loops = parse_structure(&st_str)?;

    // Resolve the preset into its two halves, then let --nn / --pk override either.
    let (preset_nn, preset_pk) = match cli.pk_params {
        PkChoice::Dp09 => (NnChoice::Mt09, DpChoice::Dp09),
        PkChoice::T04Dp03 => (NnChoice::Turner2004, DpChoice::Dp03),
        PkChoice::T04Dp09 => (NnChoice::Turner2004, DpChoice::Dp09),
        PkChoice::A07Dp03 => (NnChoice::Andronescu2007, DpChoice::Dp03),
    };
    let nn = cli.nn.unwrap_or(preset_nn);
    let pk = cli.pk.unwrap_or(preset_pk);

    let base = match nn {
        NnChoice::Mt09 => ViennaRNA::from_andrunescu_params(&RNA_MT09),
        NnChoice::Andronescu2007 => ViennaRNA::from_andrunescu_params(&RNA_ANDRONESCU_2007),
        NnChoice::Turner2004 => ViennaRNA::from_thermo_params(&RNA_TURNER_2004, cli.celsius),
        NnChoice::Extended => ViennaRNA::from_thermo_params(&RNA_EXTENDED, cli.celsius),
    };
    let model = match pk {
        DpChoice::Dp09 => base.with_pseudoknot_params(RNA_DP09),
        DpChoice::Dp03 => base.with_pseudoknot_params(RNA_DP03),
    };

    if cli.show_params {
        // DP09 = MT09 + dp09 is the only trained pairing reachable here. The
        // paper's DP03 is MT99 + dp03, and MT99 is not carried by this crate.
        let trained = matches!((nn, pk), (NnChoice::Mt09, DpChoice::Dp09));
        eprintln!(
            "nn={nn:?} pk={pk:?}{}",
            if trained {
                ""
            } else {
                "  [WARNING: these sets were not trained together]"
            }
        );
    }

    if cli.verbose {
        let mut total = 0i32;
        for (i, lp) in loops.iter().enumerate() {
            let e = model.energy_of_pseudo_loop(&seq, lp)?;
            total += e;
            println!("{:3}  {:8.4}  {lp}", i, e as f64 / 100.0);
        }
        println!("{}", "-".repeat(60));
        println!("sum  {:8.4}", total as f64 / 100.0);
    } else {
        let energy = model.energy_of_pseudoknotted_structure(&seq, &loops)?;
        println!("{:.4}", energy as f64 / 100.0);
    }
    Ok(())
}
