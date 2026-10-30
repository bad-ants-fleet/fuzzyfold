use std::sync::Arc;
use std::io::Cursor;
use std::time::Duration;
use std::hint::black_box;
use criterion::Criterion;
use criterion::SamplingMode;
use criterion::criterion_group;
use criterion::criterion_main;
use rand::SeedableRng;
use rand::rngs::StdRng;

use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_energy::Fuzzyfold;
use ff_kinetics::Arrhenius;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::ShiftMoves;
use ff_kinetics::SSA;
use ff_kinetics::MoveSet;

mod common;
use common::{data, load_inputs};

/// One registry per input structure, from the blocks of a macrostate file
/// (each block starts with a '>' line).
fn load_registries(len: usize, emodel: &Arc<Fuzzyfold>) -> Vec<MacrostateRegistry<Fuzzyfold>> {
    let path = data(&format!("benchmark_macrostates_len{len}.ms"));
    let text = std::fs::read_to_string(&path).expect("Cannot open macrostate file");
    let mut blocks: Vec<String> = Vec::new();
    for line in text.lines() {
        if line.trim_start().starts_with('>') || blocks.is_empty() {
            blocks.push(String::new());
        }
        let block = blocks.last_mut().unwrap();
        block.push_str(line);
        block.push('\n');
    }
    let inputs = load_inputs(len);
    assert_eq!(blocks.len(), inputs.len(),
        "Macrostate block count must match input trajectory count");
    blocks.iter().zip(&inputs).map(|(block, (seq, _))| {
        let mut registry = MacrostateRegistry::from((seq.clone(), emodel.clone()));
        registry.insert_from_reader(Cursor::new(block.as_bytes()), &path, 0)
            .expect("Failed to parse macrostate block");
        registry
    }).collect()
}

/// Simulations that classify the structure after every step: the cost of
/// `current_structure()` and `classify()` on top of the simulation. Every
/// iteration starts from the same seed.
fn bench_classify_dbv(c: &mut Criterion) {
    let emodel = Arc::new(Fuzzyfold::default());
    let rmodel = Arrhenius::new(emodel.temperature(), 1.0, None, None);
    let mut group = c.benchmark_group("Classify macrostates (DBV) during SSA");
    group.sampling_mode(SamplingMode::Flat);
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(10));

    for len in [100, 1000] {
        let inputs = load_inputs(len);
        let registries = load_registries(len, &emodel);
        group.bench_function(format!("dbv_len_{len:04}"), |b| {
            b.iter(|| {
                let mut rng = StdRng::seed_from_u64(42);
                let mut classifications = 0usize;
                for ((seq, pt), registry) in inputs.iter().zip(&registries) {
                    let moves = LoopNeighbors::try_from((seq.clone(), pt, emodel.clone(), ShiftMoves::NONE))
                        .expect("Failed to build LoopNeighbors");
                    let mut simulator = SSA::from((moves, rmodel));
                    simulator.simulate(&mut rng, black_box(10.0), |_t, _ti, _fl, moves| {
                        classifications += registry.classify(&moves.current_structure());
                        true
                    });
                }
                black_box(classifications)
            })
        });
    }
    group.finish();
}

/// The cost of classify() alone, without a simulation.
fn bench_classify_isolated(c: &mut Criterion) {
    let emodel = Arc::new(Fuzzyfold::default());
    let mut group = c.benchmark_group("Isolated classify() calls");

    let structures: Vec<DotBracketVec> = load_inputs(50).iter()
        .map(|(_, pt)| DotBracketVec::try_from(pt).unwrap())
        .collect();
    let registries = load_registries(50, &emodel);

    group.bench_function("isolated_dbv", |b| {
        b.iter(|| {
            let mut sum = 0usize;
            for (s, registry) in structures.iter().zip(&registries) {
                sum += registry.classify(black_box(s));
            }
            black_box(sum)
        })
    });
    group.finish();
}

criterion_group!(benches, bench_classify_dbv, bench_classify_isolated);
criterion_main!(benches);
