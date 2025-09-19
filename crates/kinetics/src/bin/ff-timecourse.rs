use clap::Args;
use clap::Parser;
use anyhow::Result;
use anyhow::bail;
use colored::*;

use std::sync::Arc;
use std::path::PathBuf;
use rayon::prelude::*;

use rand::rng;
use structure::PairTable;
use structure::DotBracketVec;
use energy::ViennaRNA;
use energy::EnergyModel;
use energy::commandline_utils::EnergyModelArguments;

use kinetics::Metropolis;
use kinetics::LoopStructure;
use kinetics::LoopStructureSSA;
use kinetics::commandline_utils::read_fasta_like_input;
use kinetics::plotting::plot_occupancy_over_time;
use kinetics::timeline::Timeline;
use kinetics::timeline::load_macrostates;

#[derive(Debug, Args)]
pub struct KineticModelParams {
    /// Metropolis rate constant (must be > 0).
    #[arg(long, default_value_t = 1e6)]
    pub k0: f64,
}

#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// The last time point of the linear scale.
    #[arg(long, default_value_t = 1e-5)]
    pub t_ext: f64,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Number of time points on the linear scale: [0..t-ext]
    #[arg(long, default_value_t = 1)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale: [t-ext..t-end]
    #[arg(long, default_value_t = 20)]
    pub t_log: usize,
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        if self.t_end <= self.t_ext {
            bail!("t_end ({}) must be greater than t_ext ({})", self.t_end, self.t_ext);
        }
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        Ok(())
    }

    pub fn get_output_times(&self) -> Vec<f64> {
        let t_end = self.t_end;
        let t_ext = self.t_ext;
        let t_lin = self.t_lin;
        let t_log = self.t_log;
        let mut times = vec![0.0];

        // Linear segments: append `t_lin` evenly spaced points
        let start = *times.last().unwrap();
        let step = t_ext / t_lin as f64;
        for i in 1..=t_lin {
            times.push(start + i as f64 * step);
        }

        // Logarithmic tail
        let start = *times.last().unwrap();
        let log_start = start.ln();
        let log_end = t_end.ln();
        for i in 1..t_log {
            let frac = i as f64 / t_log as f64;
            let value = (log_start + frac * (log_end - log_start)).exp();
            times.push(value);
        }
        times.push(t_end);

        times
    }
}

#[derive(Debug, Parser)]
#[command(name = "ff-simulate")]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

    #[arg(short, long, default_value_t = 1)]
    num_sims: usize,

    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    #[command(flatten, next_help_heading = "Simulation parameters")]
    simulation: TimelineParameters,

    #[command(flatten, next_help_heading = "Kinetic model parameters")]
    kinetics: KineticModelParams,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    cli.simulation.validate()?;

    // --- Build simulator ---
    let emodel = ViennaRNA::default();
    let rmodel = Metropolis::new(emodel.temperature(), cli.kinetics.k0);

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
    let pairings = PairTable::try_from(&structure)?;

    if let Some(h) = header {
        println!("{}", h.yellow())
    }
    println!("{}", sequence);

    println!("Output after {} simulations: \n - {:?}\n - {:?}",
        cli.num_sims, cli.kinetics, cli.simulation);

    let times = cli.simulation.get_output_times();
    let registry = load_macrostates(
        &cli.macrostates,
        &sequence, 
        Some(&emodel));

    println!("Macrostates:\n{}", registry.iter()
        .map(|(_, m)| format!(" - {} {:6}", m.name(), 
                if let Some(e) = m.energy() {
                    format!("{:6.2}", e)
                } else { "".to_string() }
        ))
        .collect::<Vec<_>>().join("\n"));

    let shared_registy = Arc::new(registry);

    let timelines: Vec<Timeline> = (0..cli.num_sims)
        .into_par_iter()
        .map(|_| {
            let registry = Arc::clone(&shared_registy);
            let mut timeline = Timeline::new(&times, registry);

            let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel)).unwrap();
            let mut simulator = LoopStructureSSA::from((loops, &rmodel));
            let mut t_idx = 0;
            simulator.simulate(
                &mut rng(), 
                cli.simulation.t_end,
                |t, tinc, _, ls| {
                    while t_idx < times.len() && t+tinc >= times[t_idx] {
                        let structure = DotBracketVec::from(ls);
                        timeline.assign_structure(t_idx, &structure);
                        t_idx += 1;
                    }
                }
            );
            timeline
        })
        .collect();
    
    // Master timeline
    let mut master = Timeline::new(&times, Arc::clone(&shared_registy));
    for timeline in timelines {
        master.merge(timeline);
    }

    println!("Calculated Timeline:\n{}", master);
    plot_occupancy_over_time(&master, "myplot.svg", cli.simulation.t_ext, cli.simulation.t_end);

    Ok(())
}



