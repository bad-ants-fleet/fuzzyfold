use rand::rng;
use clap::Parser;
use anyhow::Result;
use anyhow::bail;
use colored::*;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Metropolis;
use ff_kinetics::Kawasaki;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy;
use ff_kinetics::SSA;

use fuzzyfold::input_parsers::read_cotr_input;
use fuzzyfold::input_parsers::read_eval_input;
use fuzzyfold::energy_parsers::EnergyModelArguments;
use fuzzyfold::kinetics_parsers::RateModelArguments;
use fuzzyfold::kinetics_parsers::RateModelKind;
//TODO: support seeded rng.

#[derive(Debug, Parser)]
#[command(version, about = "Stochastic Simulation Algorithm for RNA folding")]
pub struct Cli {
    /// Input file (FASTA-like), or "-" for stdin
    #[arg(value_name = "INPUT", default_value = "-")]
    input: String,

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
    let cli = Cli::parse();

    // --- Build simulator ---
    let emodel = cli.energy.build_model();


    let (header, sequence, structure) =
        if cli.t_ext.is_some() {
            read_cotr_input(&cli.input)?
        } else {
            match read_eval_input(&cli.input) {
                Ok(v) => v,
                Err(e) => return Err(anyhow::anyhow!("{e} (or use --t-ext?)")),
            }
        };

    let pairings = PairTable::try_from(&structure)?;
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

    let times: Vec<f64> = if let Some(t_ext) = cli.t_ext {
        let mut v = vec![t_ext; sequence.len() - structure.len()];
        v.push(cli.t_end);
        v
    } else { 
        vec![cli.t_end] 
    };

    let clk = cli.kinetics;
    match (clk.rate_model, clk.k3ws.is_some(), clk.k4ws.is_some()) {
        (RateModelKind::Metropolis, false, false) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, shift_policy::NoShift))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, &times, cli.silent, cli.num_steps);
        },
        (RateModelKind::Metropolis, true, false) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, shift_policy::ThreeWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, &times, cli.silent, cli.num_steps);
        },
        (RateModelKind::Metropolis, false, true) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, shift_policy::FourWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, &times, cli.silent, cli.num_steps);
        },
        (RateModelKind::Metropolis, true, true) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, shift_policy::ThreeAndFour))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, &times, cli.silent, cli.num_steps);
        },
        (RateModelKind::Kawasaki, false, false) => {
            let rmodel = Kawasaki::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, shift_policy::NoShift))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, &times, cli.silent, cli.num_steps);
        },
        (RateModelKind::Kawasaki, _, _) => {
            bail!("Shift moves are only available for the Metropolis model.")
        },
    }

    Ok(())
}


fn run_simulator<'a, W: Walker + std::fmt::Display, K: RateModel>(
    moves: W,
    rmodel: &'a K,
    times: &[f64],
    silent: bool,
    num_steps: usize,
)
where
    SSA<'a, W, K>: From<(W, &'a K)>,
{
    let mut simulator = SSA::from((moves, rmodel));

    match (silent, num_steps) {
        (true, 0) => {
            simulator.co_simulate(&mut rng(), times, |_, _, _, _|  true);
            let cs = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", cs, en as f64 / 100.);
        },
        (false, 0) => {
            simulator.co_simulate(&mut rng(), times, 
                |t, tinc, flux, w| {
                    println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                        w,
                        w.current_energy() as f64 / 100.,
                        t,
                        tinc,
                        1.0 / flux,
                    );
                    true
                });
        },
        (true, _) => {
            if times.len() > 1 {
                panic!("no support for co-transcriptional num-steps.")
            }
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, |_, _, _, _| {
                steps += 1;
                steps < num_steps 
            });
            let ls = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", ls, en as f64 / 100.);
        },
        (false, _) => {
            if times.len() > 1 {
                panic!("no support for co-transcriptional num-steps.")
            }
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, 
                |t, tinc, flux, w| {
                    println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                        w,
                        w.current_energy() as f64 / 100.,
                        t,
                        tinc,
                        1.0 / flux,
                    );
                    steps += 1;
                    steps < num_steps 
                });
        },
    }
}
