use std::sync::Arc;
use std::hint::black_box;
use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;

use ff_energy::Fuzzyfold;
use ff_kinetics::shift_policy::*;
use ff_kinetics::LoopNeighbors;

mod common;
use common::load_inputs;

const LENGTHS: [usize; 2] = [50, 100];

pub fn gen_neighbors(c: &mut Criterion) {
    let emodel = Arc::new(Fuzzyfold::default());
    let mut group = c.benchmark_group("Generate neighborhoods");
    group.measurement_time(std::time::Duration::from_secs(20));
    group.sample_size(10);

    for len in LENGTHS {
        let inputs = load_inputs(len);
        group.bench_function(format!("enumerate_len_{len:04}"), |b| {
            b.iter(|| {
                let mut count = 0usize;
                for (seq, pt) in &inputs {
                    let mut moves = LoopNeighbors::try_from((seq.clone(), pt, emodel.clone(), NoShift))
                        .expect("Failed to build LoopNeighbors");
                    moves.generate_neighbors(0, 5, |_, _| {
                        count += 1;
                    });
                }
                black_box(count)
            })
        });
    }
    group.finish();
}

criterion_group!(benches, gen_neighbors);
criterion_main!(benches);
