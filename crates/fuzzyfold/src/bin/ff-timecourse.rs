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
use ff_kinetics::AddDelMoves;
use ff_kinetics::AddDelShiftMoves;
use ff_kinetics::Walker;
use ff_kinetics::SSA;
use ff_kinetics::timeline::Timeline;
use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::RateModel;

use fuzzyfold::kinetics_parsers::RateModelKind;
use fuzzyfold::input_parsers::read_fasta_like_input;
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

    /// Use shift moves.
    #[arg(long)]
    pub shift_moves: bool, 

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    /// Backup/Store timeline in this file.
    #[arg(long, value_name = "FILE")]
    timeline: Option<PathBuf>,

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
    let emodel = cli.energy.build_model();

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;

    let name = if let Some(h) = header {
        println!("{}", h.yellow());
        h.strip_prefix('>')
            .and_then(|s| s.split_whitespace().next())
            .unwrap_or("anonymous")
            .to_string()
    } else {
        println!("{}", ">anonymous".yellow());
        "anonymous".to_string()
    };
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation, cli.energy);

    let times = cli.simulation.get_output_times();
    let mut registry = MacrostateRegistry::from((&sequence, &emodel));
    registry.insert_files(&cli.macrostates)?;

    println!("Macrostates:\n{}", registry.iter()
        .map(|(_, m)| format!(" - {} {:6.2}", m.name(), m.ensemble_energy().unwrap_or(0.0)))
        .collect::<Vec<_>>().join("\n"));

    let shared_registry = Arc::new(registry);

    // If timeline.json exists, reload instead of starting empty
    let mut master = if let Some(path) = &cli.timeline {
        if Path::new(path).exists() {
            println!("Loading existing timeline from: {}", path.display());
            Timeline::from_file(path, &times, Arc::clone(&shared_registry))?
        } else {
            println!("A new timeline file will be created: {}", 
                path.display());
            Timeline::new(&times, Arc::clone(&shared_registry))
        }
    } else {
        Timeline::new(&times, Arc::clone(&shared_registry))
    };

    match (cli.shift_moves, cli.kinetics.model) {
        (true, RateModelKind::Metropolis) => {
            let rmodel = cli.kinetics.build_metropolis_model(emodel.temperature());
            let moves = AddDelShiftMoves::try_from((&sequence, &pairings, &emodel))
                .map_err(|e| anyhow::anyhow!(
                        "failed to construct AddDelMoves: {:?}", e
                ))?;
            let timelines: Vec<_> = run_timecourse::<AddDelShiftMoves<_>, _, _>(
                moves,
                &times,
                &rmodel,
                &cli,
                Arc::clone(&shared_registry),
            ).collect();

            for timeline in timelines {
                master.merge(timeline);
            }
            finalize(master, &cli, name);
        },
        (false, RateModelKind::Metropolis) => {
            let rmodel = cli.kinetics.build_metropolis_model(emodel.temperature());
            let moves = AddDelMoves::try_from((&sequence, &pairings, &emodel))
                .map_err(|e| anyhow::anyhow!(
                        "failed to construct AddDelMoves: {:?}", e
                ))?;
            let timelines: Vec<_> = run_timecourse::<AddDelMoves<_>, _, _>(
                moves,
                &times,
                &rmodel,
                &cli,
                Arc::clone(&shared_registry),
            ).collect();

            for timeline in timelines {
                master.merge(timeline);
            }
            finalize(master, &cli, name);
        },
        (true, RateModelKind::Kawasaki) => {
            let rmodel = cli.kinetics.build_kawasaki_model(emodel.temperature());
            let moves = AddDelShiftMoves::try_from((&sequence, &pairings, &emodel))
                .map_err(|e| anyhow::anyhow!(
                        "failed to construct AddDelMoves: {:?}", e
                ))?;
            let timelines: Vec<_> = run_timecourse::<AddDelShiftMoves<_>, _, _>(
                moves,
                &times,
                &rmodel,
                &cli,
                Arc::clone(&shared_registry),
            ).collect();

            for timeline in timelines {
                master.merge(timeline);
            }
            finalize(master, &cli, name);
        },
        (false, RateModelKind::Kawasaki) => {
            let rmodel = cli.kinetics.build_kawasaki_model(emodel.temperature());
            let moves = AddDelMoves::try_from((&sequence, &pairings, &emodel))
                .map_err(|e| anyhow::anyhow!(
                        "failed to construct AddDelMoves: {:?}", e
                ))?;
            let timelines: Vec<_> = run_timecourse::<AddDelMoves<_>, _, _>(
                moves,
                &times,
                &rmodel,
                &cli,
                Arc::clone(&shared_registry),
            ).collect();

            for timeline in timelines {
                master.merge(timeline);
            }
            finalize(master, &cli, name);
        },
    }
    Ok(())
}

fn finalize<'a, E: EnergyModel>(master: Timeline<'a, E>, cli: &Cli, name: String) {
    println!("Final Timeline:\n{}", master);
    plot_occupancy_over_time(&master, &format!("ff_{}.svg", name), cli.simulation.t_ext, cli.simulation.t_end);

    if let Some(path) = &cli.timeline {
        let serial = master.to_serializable();
        let json = to_string_pretty(&serial).unwrap();
        fs::write(path, json).unwrap();
    }
}

fn run_timecourse<'t, W, K, E>(
    moves: W,
    times: &'t [f64],
    rmodel: &'t K,
    cli: &Cli,
    registry: Arc<MacrostateRegistry<'t, E>>,
) -> impl ParallelIterator<Item = Timeline<'t, E>>
where
    W: Walker + Clone + Send + Sync,
    K: RateModel + Sync,
    E: EnergyModel + Sync,
    SSA<'t, W, K>: From<(W, &'t K)>,
{
    println!("Simulation progress:");
    let pb = ProgressBar::new(cli.num_sims as u64);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    (0..cli.num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _| {
                let registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(times, registry);

                let mut simulator = SSA::from((moves.clone(), rmodel));
                let mut t_idx = 0;
                simulator.simulate(
                    &mut rng(),
                    cli.simulation.t_end,
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


