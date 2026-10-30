use rustc_hash::FxHashMap;
use std::sync::Arc;
use std::path::PathBuf;
use std::io::{BufWriter, Write};
use std::fs::File;

use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use ff_structure::DotBracketVec;
use ff_structure::MultiPairTable;
use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_kinetics::Elongation;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::SSA;
use ff_kinetics::RateModel; 

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimeArguments;
use ff_kinetics::OutputTimes;

#[derive(Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimeArguments,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let emodel = Arc::new(cli.energy.build_model()?);
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) =
        if cli.simulation.t_ext.is_some() {
            read_cotr_input(&cli.input, is_rna)?
        } else {
            match read_eval_input(&cli.input, is_rna) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };
    cli.energy.warn_circular(&sequence);
    let num_ext = sequence.len() - structure.len();
    let pairings = MultiPairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow());
    }
    println!("{}", sequence);
 
    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let output = cli.simulation.output_times(num_ext)?;
    let schedule = cli.simulation.schedule(num_ext, &output)?;
    let times = output.times();
    
    let edits = rmodel.move_types().a_to_i;
    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.move_types()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    let trajectories: Vec<Vec<State>> =
        run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64, &output).collect();

    let mut master: Vec<FxHashMap<State, usize>> = vec![FxHashMap::default(); times.len()];
    for trajectory in trajectories {
        for (i, state) in trajectory.into_iter().enumerate() {
             *master[i].entry(state).or_insert(0) += 1;
        }
    }
    println!("{}", "Finished simulations!".red());

    let drf_path = cli.output.with_extension("drf");
    let idict = build_idict(&master, sequence.len());
    let mut writer = BufWriter::new(File::create(drf_path.clone())?);
    // With A-to-I edits, the sequence of each structure is a column.
    if edits {
        writeln!(writer, "id time occupancy structure energy sequence")?;
    } else {
        writeln!(writer, "id time occupancy structure energy")?;
    }

    for (t_idx, ensemble) in master.iter().enumerate() {
        let t = times[t_idx];
        let mut entries: Vec<(&State, usize)> = ensemble.iter().map(|(s, &c)| (s, c)).collect();
        entries.sort_by_key(|&((_, _, e), _)| *e);

        for (state @ (seq, s, en), occu) in entries {
            let id = idict[&id_key(state, sequence.len())];
            write!(
                writer,
                "{:5} {:03.3} {:5} {} {:6.2}",
                id, t, occu as f64 / cli.num_sims as f64, s, *en as f64 / 100.0
            )?;
            if edits {
                write!(writer, " {}", &seq.to_string()[..s.len()])?;
            }
            writeln!(writer)?;
        }
    }
    println!("Wrote drf file: {}", format!("{}", drf_path.display()).green()); 

    Ok(())
}

/// A structure with its sequence (full length, with edits) and energy.
type State = (NucleotideVec, DotBracketVec, i32);

/// One id per sequence and structure, the structure padded to full length:
/// a structure keeps its id when the strand grows, not when it is edited.
fn build_idict(master: &[FxHashMap<State, usize>], seqlen: usize) -> FxHashMap<String, usize> {
    let mut idict: FxHashMap<String, usize> = FxHashMap::default();
    for states in master {
        for state in states.keys() {
            let next_id = idict.len();
            idict.entry(id_key(state, seqlen)).or_insert(next_id);
        }
    }
    idict
}

fn id_key((seq, s, _): &State, seqlen: usize) -> String {
    let mut key = s.to_string();
    key.extend(std::iter::repeat_n('.', seqlen - s.len()));
    key.push(' ');
    key.push_str(&seq.to_string());
    key
}

fn run_timecourse<W, K>(
    moves: W,
    rmodel: K,
    schedule: &[f64],
    num_sims: u64,
    output: &OutputTimes,
) -> impl ParallelIterator<Item = Vec<State>>
where
    W: Elongation + Clone + Send + Sync,
    K: RateModel + Clone + Send + Sync,
{
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _| {
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let trajectory = simulator.sample(&mut rng(), schedule, output,
                    |w| (w.current_sequence().clone(), w.current_structure(), w.current_energy()));
                pb.inc(1);
                trajectory
            },
        )
}
