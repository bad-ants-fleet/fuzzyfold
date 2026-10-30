use clap::Parser;
use anyhow::Result;
use colored::*;
use rand::rng;
use std::sync::Arc;

use fuzzyfold::structure::MultiPairTable;
use fuzzyfold::energy::EnergyModel;
use fuzzyfold::kinetics::RateModel;
use fuzzyfold::kinetics::Elongation;
use fuzzyfold::kinetics::LoopNeighbors;
use fuzzyfold::kinetics::SSA;
use fuzzyfold::kinetics::extension_schedule;

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::input_parsers::read_fasta_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
//TODO: support seeded rng.

#[derive(Debug, Parser)]
#[command(version, about = "One stochastic folding trajectory, also co-transcriptional (--t-ext).", after_help = fuzzyfold::CITATION)]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    /// Read input as standard FASTA (DNA alphabet, multi-line sequence, no structure)
    #[arg(long)]
    fasta: bool,

    /// Simulation time at the full-length molecule.
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

    /// Optional simulation time per nucleotide during co-transcriptional folding.
    #[arg(long)]
    t_ext: Option<f64>,

    /// Set the number of simulation steps (ignore t_end!!)
    #[arg(long, hide = true, default_value_t = 0)]
    num_steps: usize,

    /// Do not print trajectory, only last structure.
    #[arg(long)]
    silent: bool, 

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    fuzzyfold::restore_sigpipe();
    let cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model()?);
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) =
        if cli.fasta {
            read_fasta_input(&cli.input, true)?  // is_rna is ignored, assuming DNA alphabet
        } else if cli.t_ext.is_some() {
            read_cotr_input(&cli.input, is_rna)?
        } else {
            match read_eval_input(&cli.input, is_rna) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };
    cli.energy.warn_circular(&sequence);

    let pairings = MultiPairTable::try_from(&structure)?;
    if let Some(h) = header {
        println!("{}", h.yellow())
    }
    println!("{} {:>8} {:>14} {:>14} {:>15}",
        sequence,
        "energy".green(),
        "arrival-time".cyan(),
        "waiting-time".cyan(),
        "mean-waiting".cyan(),
    );

    let num_ext = sequence.len() - structure.len();
    let end = cli.t_ext.map_or(0.0, |t| t * num_ext as f64) + cli.t_end;
    let times = extension_schedule(cli.t_ext, num_ext, end).map_err(anyhow::Error::msg)?;

    let moves = LoopNeighbors::try_from((sequence, &pairings, emodel, rmodel.move_types()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    run_simulator(moves, rmodel, &times, cli.silent, cli.num_steps);

    Ok(())
}


fn run_simulator<W: Elongation + std::fmt::Display, K: RateModel>(
    moves: W,
    rmodel: K,
    times: &[f64],
    silent: bool,
    num_steps: usize,
)
{
    let mut simulator = SSA::from((moves, rmodel));
    // When an edit changes the sequence, it is printed (as far as it is
    // transcribed) before the structure, with the same columns.
    let mut sequence = simulator.move_set().current_sequence().clone();
    let mut edited = move |w: &W| {
        let changed = w.current_sequence() != &sequence;
        if changed {
            sequence = w.current_sequence().clone();
        }
        changed.then(|| sequence.to_string()[..w.current_structure_length()].to_string())
    };
    let mut print_step = |t: f64, tinc: f64, flux: f64, w: &W| {
        let columns = format!("{:8.2} {:14.8e} {:14.8e} {:15.8e}",
            w.current_energy() as f64 / 100., t, tinc, 1.0 / flux);
        if let Some(seq) = edited(w) {
            println!("{} {}", seq.yellow(), columns);
        }
        println!("{} {}", w, columns);
    };

    match (silent, num_steps) {
        (true, 0) => {
            simulator.co_simulate(&mut rng(), times, |_, _, _, _|  true);
        },
        (false, 0) => {
            simulator.co_simulate(&mut rng(), times, |t, tinc, flux, w| {
                print_step(t, tinc, flux, w);
                true
            });
        },
        (silent, _) => {
            if times.len() > 1 {
                panic!("no support for co-transcriptional num-steps.")
            }
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, |t, tinc, flux, w| {
                if !silent {
                    print_step(t, tinc, flux, w);
                }
                steps += 1;
                steps < num_steps
            });
        },
    }
    if silent {
        let w = simulator.move_set();
        let energy = format!("{:8.2}", w.current_energy() as f64 / 100.);
        if let Some(seq) = edited(w) {
            println!("{} {}", seq.yellow(), energy);
        }
        println!("{} {}", simulator.current_structure(), energy);
    }
}
