use std::fs::File;
use std::hint::black_box;
use std::io::BufRead;
use std::io::BufReader;
use criterion::Criterion;
use criterion::BatchSize;
use criterion::criterion_group;
use criterion::criterion_main;

use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA;
use ff_kinetics::AddDelMoves;

const INPUT_L50: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len50.vrna");

const INPUT_L100: &str = concat!(env!("CARGO_MANIFEST_DIR"), 
    "/benches/data/benchmark_random_structures_len100.vrna");

struct BenchCase {
    name: &'static str,
    path: &'static str,
}

const CASES: &[BenchCase] = &[
    BenchCase { name: "len_0050", path: INPUT_L50 },
    BenchCase { name: "len_0100", path: INPUT_L100 },
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

pub fn gen_neighbors(c: &mut Criterion) {
    let emodel = ViennaRNA::default();
    let mut group = c.benchmark_group("Generate neighborhoods");
    group.measurement_time(std::time::Duration::from_secs(20));
    group.sample_size(10);

    for case in CASES {
        let inputs = load_raw_inputs(case.path); 
        group.bench_function(format!("enumerate_{}", case.name), |b| {
            b.iter_batched(
                || &inputs, 
                |inputs| {
                    let mut count = 0usize;
                    for (seq, pt) in inputs {
                        let mut moves = AddDelMoves::try_from((seq, pt, &emodel))
                            .expect("failed to build loop table");
                        moves.generate_neighbors(0, 5, |_, _| {
                            count += 1;
                        });
                    }
                    black_box(count)
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

criterion_group!(benches, gen_neighbors);
criterion_main!(benches);
