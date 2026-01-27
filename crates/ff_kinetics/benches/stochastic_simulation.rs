use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::hint::black_box;
use criterion::BatchSize;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::Criterion;
use rand::SeedableRng;
use rand::rngs::StdRng;

use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;
use ff_kinetics::Metropolis;
use ff_kinetics::AddDelMoves;
use ff_kinetics::AddDelShiftMoves;
use ff_kinetics::SSA;

const INPUT_L50: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len50.vrna");

const INPUT_L100: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len100.vrna");

const INPUT_L250: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len250.vrna");

const INPUT_L500: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len500.vrna");

const INPUT_L750: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len750.vrna");

const INPUT_L1000: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len1000.vrna");

const INPUT_L2500: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len2500.vrna");

struct BenchCase {
    name: &'static str,
    path: &'static str,
}

const CASES: &[BenchCase] = &[
    BenchCase { name: "len_0050", path: INPUT_L50 },
    BenchCase { name: "len_0100", path: INPUT_L100 },
    BenchCase { name: "len_0250", path: INPUT_L250 },
    BenchCase { name: "len_0500", path: INPUT_L500 },
    BenchCase { name: "len_0750", path: INPUT_L750 },
    BenchCase { name: "len_1000", path: INPUT_L1000 },
    BenchCase { name: "len_2500", path: INPUT_L2500 },
];

fn load_raw_inputs(path: &str) -> Vec<(NucleotideVec, PairTable)> {
    let file = File::open(path).expect("Cannot open input file");
    let reader = BufReader::new(file);
    let mut inputs = Vec::new();
    let mut lines = reader.lines();
    while let Some(Ok(header)) = lines.next() {
        assert!(header.starts_with('>'), "Malformed benchmarking input.");
        let seq = lines.next().unwrap().unwrap();
        let dbr = lines.next().unwrap().unwrap();
        let seq = NucleotideVec::try_from(seq.as_str()).unwrap();
        let pt  = PairTable::try_from(dbr.as_str()).unwrap();
        inputs.push((seq, pt));
    }
    inputs
}

fn simulate_benchmark(c: &mut Criterion) {
    let emodel = ViennaRNA::default();
    let rmodel = Metropolis::new(emodel.temperature(), 1.0, 1.0);
    let mut group = c.benchmark_group("Seeded stochastic simulations.");
    group.measurement_time(std::time::Duration::from_secs(50)); 

    for case in CASES {
        let inputs = load_raw_inputs(case.path); 
        let mut rng = StdRng::seed_from_u64(42);
        group.bench_function(format!("simulate_{}", case.name), |b| {
            b.iter_batched(
                || &inputs, 
                |inputs| {
                    for (seq, pt) in inputs {
                        let moves = AddDelMoves::try_from((seq, pt, &emodel))
                            .expect("failed to build loop table");
                        let mut simulator = SSA::from((moves, &rmodel));

                        simulator.simulate(
                            &mut rng, 
                            black_box(10.0), 
                            |_t, _ti, _fl, _w| { true }
                        );
                    }
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

fn simulate_shift_benchmark(c: &mut Criterion) {
    let emodel = ViennaRNA::default();
    let rmodel = Metropolis::new(emodel.temperature(), 1.0, 1.0);
    let mut group = c.benchmark_group("Seeded stochastic simulations with shift moves.");
    group.measurement_time(std::time::Duration::from_secs(50)); 

    for case in CASES {
        let inputs = load_raw_inputs(case.path); 
        let mut rng = StdRng::seed_from_u64(42);
        group.bench_function(format!("simulate_shift_{}", case.name), |b| {
            b.iter_batched(
                || &inputs, 
                |inputs| {
                    for (seq, pt) in inputs {
                        let moves = AddDelShiftMoves::try_from((seq, pt, &emodel))
                            .expect("failed to build loop table");
                        let mut simulator = SSA::from((moves, &rmodel));

                        simulator.simulate(
                            &mut rng, 
                            black_box(10.0), 
                            |_t, _ti, _fl, _w| { true }
                        );
                    }
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

criterion_group!(benches, simulate_benchmark, simulate_shift_benchmark);
criterion_main!(benches);

