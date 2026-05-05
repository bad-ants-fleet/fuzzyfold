use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::BufWriter;
use std::sync::Arc;
use std::path::Path;
use std::path::PathBuf;

use rayon::prelude::*;
use rand::rng;
use colored::*;
use clap::Parser;
use anyhow::Result;
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

    #[arg(short, long)]
    title: Option<String>,

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
    // Verbose Output
    println!("{:>4} {:<10} {} {:>5} {:>8}",
        "ID",
        "Macrostate".cyan(), format!("{}", sequence).yellow(), "Size", "Energy");
    for (id, m) in macrostates.iter() {
        if m.name() == "Unassigned" {
            println!("{:4} {:<10}", 0, m.name());
            continue
        }
        println!("{:4} {:<10} {:<} {:>5} {:>8.2}",
            id, 
            m.name(),
            m.get_lowest_microstate().unwrap(),
            m.len(),
            m.ensemble_energy().unwrap());
    }
    let shared_macrostates = Arc::new(macrostates);

    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");
    let nxy_path = cli.output.with_extension("nxy");

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

    println!("{}", "Finished simulations!".red());

    // save / print / plot.
    let mut writer = BufWriter::new(File::create(nxy_path.clone())?);
    write!(writer, "{}", master)?;
    println!("Wrote nxy file: {}", format!("{}",nxy_path.display()).green());

    let serial = master.to_serializable();
    let json = to_string_pretty(&serial).unwrap();
    fs::write(tln_path.clone(), json).unwrap();
    println!("Wrote tln file: {}", tln_path.display());

    let numsim = master.points[0].counter;
    let title = cli.title.unwrap_or({
        format!("ff-timecourse ({} simulations)", 
        {
            if numsim >= 10000 {
                let s = numsim.to_string();
                let mut out = String::new();
                for (i, c) in s.chars().rev().enumerate() {
                    if i > 0 && i % 3 == 0 {
                        out.push('_');
                    }
                    out.push(c);
                }
                out.chars().rev().collect::<String>()
            } else { 
                numsim.to_string()
            }
        })
    });

    plot_occupancy_over_time(&master, svg_path.clone(), &title, cli.simulation.t_ext, cli.simulation.t_end);
    println!("Plotted svg file: {}", svg_path.display());

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
