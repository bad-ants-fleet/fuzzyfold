//! A fixed workload for instruction counts: seeded simulations on the
//! benchmark structures, without shift moves (100 nt, to time 10) or with
//! three- and four-way shifts (250 nt, 3000 steps each). Prints a sum of the
//! energies and steps (identical results are part of every comparison).
//!
//! cargo build --release -p ff_kinetics --example kinetics_workload
//! valgrind --tool=callgrind target/release/examples/kinetics_workload shifts

use std::sync::Arc;
use rand::SeedableRng;
use rand::rngs::StdRng;
use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, SSA};

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../ff_energy/benches/data");

fn run(shifts: bool, len: usize, count: usize, t_max: f64, steps: i64) -> i64 {
    let model = Arc::new(Fuzzyfold::default());
    let s = shifts.then_some(1.0);
    let rates = Arrhenius::new(model.temperature(), 1.0, s, s, None);
    let mut rng = StdRng::seed_from_u64(42);
    let input = std::fs::read_to_string(format!("{DATA}/benchmark_random_structures_len{len}.vrna")).unwrap();
    let lines: Vec<&str> = input.lines().collect();
    let mut sum = 0;
    for c in lines.chunks(3).take(count) {
        let seq = NucleotideVec::try_from(c[1]).unwrap();
        let start = MultiPairTable::try_from(c[2]).unwrap();
        let moves = LoopNeighbors::try_from((seq, &start, model.clone(), rates.move_types())).unwrap();
        let mut ssa = SSA::from((moves, rates));
        let mut n = 0i64;
        ssa.simulate(&mut rng, t_max, |_, _, _, w| { n += 1; sum += w.current_energy() as i64; n < steps });
        sum += n;
    }
    sum
}

fn main() {
    let what = std::env::args().nth(1).expect("Usage: kinetics_workload <no-shifts|shifts>");
    let sum = match what.as_str() {
        "no-shifts" => run(false, 100, 100, 10.0, i64::MAX),
        "shifts" => run(true, 250, 20, 1e6, 3000),
        _ => panic!("Usage: kinetics_workload <no-shifts|shifts>"),
    };
    println!("{what} {sum}");
}
