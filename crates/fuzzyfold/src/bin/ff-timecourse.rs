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

use ff_structure::MultiPairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Elongation;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::SSA;
use ff_kinetics::timecourse::Timecourse;
use ff_kinetics::MacrostateRegistry;

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

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[arg(short, long)]
    title: Option<String>,

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
                                                                            //
    let mut macrostates = MacrostateRegistry::from((sequence.clone(), emodel.clone()));
    macrostates.insert_files(&cli.macrostates, cli.simulation.t_ext.is_some())?;
    // Verbose Output
    println!("{:>4} {:<10} {} {:>5} {:>8}",
        "ID",
        "Macrostate".cyan(), format!("{}", sequence).yellow(), "Size", "Energy");
    for (id, (l, m)) in macrostates.iter() {
        if *l != sequence.len() {
            continue
        }
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

    // Output paths 
    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");
    let nxy_path = cli.output.with_extension("nxy");

    // If the .tln file exists, add to it instead of starting empty
    let mut master = 
        if Path::new(&tln_path).exists() {
            println!("Loading existing time course from: {}", tln_path.display());
            Timecourse::read_tln(&tln_path, times, &macrostates)?
        } else {
            println!("A new time course file will be created: {}", tln_path.display());
            Timecourse::new(times, macrostates.len())
        };

    let moves = LoopNeighbors::try_from((sequence.clone(), &pairings, emodel, rmodel.shift_moves()))
        .map_err(|e| anyhow::anyhow!("failed to build the move set: {e}"))?;
    master.merge(&run_timecourse(moves, rmodel, &schedule, cli.num_sims as u64,
        &macrostates, &output));

    println!("{}", "Finished simulations!".red());

    // save / print / plot.
    let mut writer = BufWriter::new(File::create(nxy_path.clone())?);
    master.write_nxy(&mut writer, &macrostates)?;
    writer.flush()?;
    println!("Wrote nxy file: {}", format!("{}",nxy_path.display()).green());    
    master.write_tln(&tln_path, &macrostates)?;
    println!("Wrote tln file: {}", tln_path.display());

    let numsim = master.num_trajectories(0);
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
    let t_sep = output.t_sep().expect("lin_log output times");
    master.plot_svg(&svg_path, &macrostates, &title, t_sep, output.end())
        .map_err(anyhow::Error::from_boxed)?;
    println!("Wrote svg file: {}", svg_path.display());

    Ok(())
}


fn run_timecourse<W, K, E>(
    moves: W,
    rmodel: K,
    schedule: &[f64],
    num_sims: u64,
    registry: &MacrostateRegistry<E>,
    output: &OutputTimes,
) -> Timecourse
where
    W: Elongation + Clone + Send + Sync,
    K: RateModel + Clone,
    E: EnergyModel,
{
    let pb = ProgressBar::new(num_sims);
    pb.set_style(
        ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"),
    );

    let empty = || Timecourse::new(output.times(), registry.len());
    (0..num_sims)
        .into_par_iter()
        .map_init(
            move || pb.clone(), // each thread gets a clone
            move |pb, _| {
                let mut simulator = SSA::from((moves.clone(), rmodel.clone()));
                let trajectory = simulator.sample(&mut rng(), schedule, output,
                    |w| registry.classify(&w.current_structure()));
                pb.inc(1);
                trajectory
            },
        )
        .fold(empty, |mut tc, trajectory| { tc.add_trajectory(&trajectory); tc })
        .reduce(empty, |mut a, b| { a.merge(&b); a })
}
