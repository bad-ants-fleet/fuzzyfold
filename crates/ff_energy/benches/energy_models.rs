use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::hint::black_box;
use criterion::BatchSize;
use criterion::Criterion;
use criterion::criterion_main;
use criterion::criterion_group;

use ff_structure::PairTable;
use ff_structure::MultiPairTable;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;
use ff_energy::NucleotideVec;

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

fn insert_quartile_separators(s: &str) -> String {
    let n = s.len();
    let q1 = n / 4;
    let q2 = n / 2;
    let q3 = 3 * n / 4;
    let mut out = String::with_capacity(n + 3);
    out.push_str(&s[..q1]);
    out.push('&');
    out.push_str(&s[q1..q2]);
    out.push('&');
    out.push_str(&s[q2..q3]);
    out.push('&');
    out.push_str(&s[q3..]);
    out
}

fn load_multi_inputs(path: &str) -> Vec<(NucleotideVec, MultiPairTable)> {
    let file = File::open(path).expect("Cannot open input file");
    let reader = BufReader::new(file);
    let mut inputs = Vec::new();
    let mut lines = reader.lines();
    while let Some(Ok(header)) = lines.next() {
        assert!(header.starts_with('>'), "Malformed benchmarking input.");
        let seq = insert_quartile_separators(&lines.next().unwrap().unwrap());
        let dbr = insert_quartile_separators(&lines.next().unwrap().unwrap());
        let seq = NucleotideVec::try_from(seq.as_str()).unwrap();
        let pt  = MultiPairTable::try_from(dbr.as_str()).unwrap();
        inputs.push((seq, pt));
    }
    inputs
}

fn bench_model_init(c: &mut Criterion) {
    c.bench_function("ViennaRNA::default()", |b| {
        b.iter(|| {
            let model = ViennaRNA::default();
            black_box(model);
        })
    });
}

fn bench_evaluation(c: &mut Criterion) {
    let emodel = ViennaRNA::default();
    let mut group = c.benchmark_group("Bulk energy evaluations.");
    group.measurement_time(std::time::Duration::from_secs(30));

    for case in CASES {
        let inputs = load_raw_inputs(case.path); 
        group.bench_function(format!("evaluate_{}", case.name), |b| {
            b.iter_batched(
                || &inputs, 
                |inputs| {
                    let mut sum = 0i32;
                    for (seq, pt) in inputs {
                        sum += emodel.energy_of_structure(
                            black_box(seq),
                            black_box(pt),
                        );
                    }
                    black_box(sum);
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

fn bench_multi_evaluation(c: &mut Criterion) {
    let emodel = ViennaRNA::default();
    let mut group = c.benchmark_group("Bulk multi-energy evaluations.");
    group.measurement_time(std::time::Duration::from_secs(30));

    for case in CASES {
        let inputs = load_multi_inputs(case.path); 
        group.bench_function(format!("multi-evaluate_{}", case.name), |b| {
            b.iter_batched(
                || &inputs, 
                |inputs| {
                    let mut sum = 0i32;
                    for (seq, pt) in inputs {
                        sum += emodel.energy_of_structure(
                            black_box(seq),
                            black_box(pt),
                        );
                    }
                    black_box(sum);
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

criterion_group!(
    benches,
    bench_model_init,
    bench_evaluation,
    bench_multi_evaluation,
);
criterion_main!(benches);

