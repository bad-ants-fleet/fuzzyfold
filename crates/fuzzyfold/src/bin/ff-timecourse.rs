use std::fs;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;
use rayon::prelude::*;
use rand::rng;
use clap::Parser;
use anyhow::Result;
use colored::*;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use serde_json::to_string_pretty;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;
use ff_kinetics::timeline::Timeline;
use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistry;

use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::TimelineParameters;

#[derive(Debug, Parser)]
#[command(version, about = "Stochastically simulated nucleic acid ensembles over time.")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: RateModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = Arc::new(cli.energy.build_model());
    let rmodel = cli.kinetics.build_model(emodel.temperature());

    let is_rna = cli.energy.dna.is_none();
    let (header, sequence, structure) = read_eval_input(&cli.input, is_rna)?;
    let pairings = PairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow());
    } 
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times();
    let mut macrostates = MacrostateRegistry::from((sequence.clone(), emodel.clone()));
    macrostates.insert_files(&cli.macrostates)?;

    println!("Macrostates:\n{}", macrostates.iter()
        .map(|(_, m)| format!(" - {} {:6.2}", m.name(), m.ensemble_energy().unwrap_or(0.0)))
        .collect::<Vec<_>>().join("\n"));

    let shared_macrostates = Arc::new(macrostates);

    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");

    // If timeline.json exists, reload instead of starting empty
    let mut master = 
        if Path::new(&tln_path).exists() {
            println!("Loading existing timeline from: {}", tln_path.display());
            Timeline::from_file(&tln_path, &times, Arc::clone(&shared_macrostates))?
        } else {
            println!("A new timeline file will be created: {}", tln_path.display());
            Timeline::new(&times, Arc::clone(&shared_macrostates))
        };

    let timelines: Vec<_> =
        match (rmodel.k3ws().is_some(), rmodel.k4ws().is_some()) {
            (false, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, NoShift))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, false) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (false, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, FourWayOnly))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
            (true, true) => {
                let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, ThreeAndFour))
                    .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
                run_timecourse(moves, rmodel, cli.simulation.t_end, cli.num_sims as u64,
                    Arc::clone(&shared_macrostates), &times).collect()
            },
        };

    for timeline in timelines {
        master.merge(timeline);
    }

    println!("Final Timeline:\n{}", master);
    plot_occupancy_over_time(&master, svg_path, cli.simulation.t_ext, cli.simulation.t_end);
    let serial = master.to_serializable();
    let json = to_string_pretty(&serial).unwrap();
    fs::write(tln_path, json).unwrap();

    Ok(())
}


fn run_timecourse<W, K, E>(
    moves: W,
    rmodel: K,
    t_end: f64,
    num_sims: u64,
    registry: Arc<MacrostateRegistry<E>>,
    times: &[f64],
) -> impl ParallelIterator<Item = Timeline<E>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Clone,
    E: EnergyModel,
    SSA<W, K>: From<(W, K)>,
{
    println!("Simulation progress:");
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
                let registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(times, registry);

                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let mut t_idx = 0;
                simulator.simulate(
                    &mut rng(),
                    t_end,
                    |t, tinc, _, w| {
                        while t_idx < times.len() && t + tinc >= times[t_idx] {
                            let structure = w.current_structure();
                            timeline.assign_structure(t_idx, &structure);
                            t_idx += 1;
                        }
                        true
                    },
                );

                pb.inc(1);
                timeline
            },
        )
}


