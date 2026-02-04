use rand::rng;
use clap::Parser;
use anyhow::Result;
use colored::*;

use ff_structure::PairTable;
use ff_energy::EnergyModel;
use ff_kinetics::RateModel;
use ff_kinetics::Metropolis;
use ff_kinetics::Kawasaki;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::*;
use ff_kinetics::SSA;

use fuzzyfold::input_parsers::read_fasta_like_input;
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

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    t_end: f64,

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

    let (header, sequence, structure) = read_fasta_like_input(&cli.input)?;
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

    let clk = cli.kinetics;
    match (clk.rate_model, clk.k3ws.is_some(), clk.k4ws.is_some()) {
        (RateModelKind::Metropolis, false, false) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, NoShift))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Metropolis, true, false) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, ThreeWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Metropolis, false, true) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, FourWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Metropolis, true, true) => {
            let rmodel = Metropolis::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, ThreeAndFour))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Kawasaki, false, false) => {
            let rmodel = Kawasaki::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, NoShift))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Kawasaki, true, false) => {
            let rmodel = Kawasaki::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, ThreeWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Kawasaki, false, true) => {
            let rmodel = Kawasaki::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, FourWayOnly))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
        (RateModelKind::Kawasaki, true, true) => {
            let rmodel = Kawasaki::new(emodel.temperature(), clk.k0, clk.k3ws, clk.k4ws);
            let moves = LoopNeighbors::try_from((&sequence, &pairings, &emodel, ThreeAndFour))
                .map_err(|e| anyhow::anyhow!("failed to construct AddDelMoves: {:?}", e))?;
            run_simulator(moves, &rmodel, cli.silent, cli.t_end, cli.num_steps);
        },
    }

    Ok(())
}


fn run_simulator<'a, W: Walker + std::fmt::Display, K: RateModel>(
    moves: W,
    rmodel: &'a K,
    silent: bool,
    t_end: f64,
    num_steps: usize,
)
where
    SSA<'a, W, K>: From<(W, &'a K)>,
{
    let mut simulator = SSA::from((moves, rmodel));
    match (silent, num_steps) {
        (true, 0) => {
            simulator.simulate(&mut rng(), t_end, |_, _, _, _|  true);
            let ls = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", ls, en as f64 / 100.);
        },
        (true, _) => {
            let mut steps = 0;
            simulator.simulate(&mut rng(), f64::MAX, |_, _, _, _| {
                steps += 1;
                steps < num_steps 
            });
            let ls = simulator.current_structure();
            let en = simulator.current_energy();
            println!("{} {:8.2}", ls, en as f64 / 100.);
        },
        (false, 0) => {
            simulator.simulate(&mut rng(), t_end, 
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
        (false, _) => {
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
