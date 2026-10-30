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
    
    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.shift_moves()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    let trajectories: Vec<Vec<(DotBracketVec, i32)>> =
        run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64, &output).collect();


    let mut master: Vec<FxHashMap<(DotBracketVec, i32), usize>> = vec![FxHashMap::default(); times.len()];    
    for trajectory in trajectories {
        for (i, db) in trajectory.into_iter().enumerate() {           
             *master[i].entry(db.clone()).or_insert(0) += 1; 
        } 
    } 
    println!("{}", "Finished simulations!".red());

    let drf_path = cli.output.with_extension("drf");
    let idict = build_idict(&master, sequence.len());
    let mut writer = BufWriter::new(File::create(drf_path.clone())?);
    //writeln!(writer, "# sequence={}", sequence)?;
    //writeln!(writer, "# num_sims={}", cli.num_sims)?;
    //writeln!(writer, "# t_ext={:?}",  cli.simulation.t_ext)?;
    //writeln!(writer, "# t_end={}", cli.simulation.t_end)?;
    writeln!(writer, "id time occupancy structure energy")?;

    for (t_idx, ensemble) in master.iter().enumerate() {
        let t = times[t_idx];

        // sorted(edict[t].items(), key=lambda x: x[1])
        let mut entries: Vec<(&DotBracketVec, i32, usize)> = ensemble
            .iter()
            .map(|((s, e), &c)| (s, *e, c))
            .collect();
        entries.sort_by_key(|&(_, e, _)| e);

        for (s, en, occu) in entries {
            let id = idict[&padded_key(s, sequence.len())];
            writeln!(
                writer,
                "{:5} {:03.3} {:5} {} {:6.2}",
                id, t, occu as f64 / cli.num_sims as f64, s, en as f64 / 100.0
            )?;
        }
    }
    println!("Wrote drf file: {}", format!("{}", drf_path.display()).green()); 

    Ok(())
}

fn build_idict(
    master: &[FxHashMap<(DotBracketVec, i32), usize>],
    seqlen: usize,
) -> FxHashMap<String, usize> {
    let mut idict: FxHashMap<String, usize> = FxHashMap::default();
    let mut next_id = 0usize;
    for structures in master {
        for (s, _) in structures.keys() {
            let key = padded_key(s, seqlen);
            idict.entry(key).or_insert_with(|| {
                let id = next_id;
                next_id += 1;
                id
            });
        }
    }
    idict
}

fn padded_key(s: &DotBracketVec, seqlen: usize) -> String {
    let s = s.to_string();
    let pad = seqlen - s.len();
    let mut padded = s;
    padded.extend(std::iter::repeat_n('.', pad));
    padded
}

fn run_timecourse<W, K>(
    moves: W,
    rmodel: K,
    schedule: &[f64],
    num_sims: u64,
    output: &OutputTimes,
) -> impl ParallelIterator<Item = Vec<(DotBracketVec, i32)>>
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
                    |w| (w.current_structure(), w.current_energy()));
                pb.inc(1);
                trajectory
            },
        )
}
