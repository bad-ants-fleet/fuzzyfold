use std::sync::Arc;
use std::time::Duration;
use std::hint::black_box;
use criterion::Criterion;
use criterion::SamplingMode;
use criterion::criterion_group;
use criterion::criterion_main;
use rand::SeedableRng;
use rand::rngs::StdRng;

use ff_energy::EnergyModel;
use ff_energy::Fuzzyfold;
use ff_kinetics::Arrhenius;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::MoveTypes;
use ff_kinetics::SSA;

mod common;
use common::load_inputs;

const LENGTHS: [usize; 3] = [100, 1000, 2500];

/// Every iteration starts from the same seed, so it simulates the same
/// trajectories.
fn simulate_benchmark(c: &mut Criterion) {
    let emodel = Arc::new(Fuzzyfold::default());
    let rmodel = Arrhenius::new(emodel.temperature(), 1.0);
    let mut group = c.benchmark_group("Seeded stochastic simulations.");
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    for len in LENGTHS {
        let inputs = load_inputs(len);
        group.bench_function(format!("simulate_len_{len:04}"), |b| {
            b.iter(|| {
                let mut rng = StdRng::seed_from_u64(42);
                for (seq, pt) in &inputs {
                    let moves = LoopNeighbors::try_from((seq.clone(), pt, emodel.clone(), MoveTypes::NONE))
                        .expect("Failed to build LoopNeighbors");
                    let mut simulator = SSA::from((moves, rmodel));
                    simulator.simulate(&mut rng, black_box(10.0), |_t, _ti, _fl, _w| true);
                }
            })
        });
    }
    group.finish();
}

fn simulate_three_way_shift_benchmark(c: &mut Criterion) {
    let emodel = Arc::new(Fuzzyfold::default());
    let rmodel = Arrhenius::new(emodel.temperature(), 1.0).with_three_way(1.0);
    let mut group = c.benchmark_group("Seeded stochastic simulations with three-way shift moves.");
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(20));

    for len in LENGTHS {
        let inputs = load_inputs(len);
        group.bench_function(format!("simulate_three_way_shift_len_{len:04}"), |b| {
            b.iter(|| {
                let mut rng = StdRng::seed_from_u64(42);
                for (seq, pt) in &inputs {
                    let moves = LoopNeighbors::try_from((seq.clone(), pt, emodel.clone(), MoveTypes { three_way: true, four_way: false, a_to_i: false }))
                        .expect("Failed to build LoopNeighbors");
                    let mut simulator = SSA::from((moves, rmodel));
                    simulator.simulate(&mut rng, black_box(10.0), |_t, _ti, _fl, _w| true);
                }
            })
        });
    }
    group.finish();
}

criterion_group!(benches, simulate_benchmark, simulate_three_way_shift_benchmark);
criterion_main!(benches);
