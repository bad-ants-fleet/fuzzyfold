use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::hint::black_box;
use criterion::Criterion;
use criterion::criterion_main;
use criterion::criterion_group;

use ff_structure::PairTable;
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
    group.measurement_time(std::time::Duration::from_secs(50)); // increase from default 5s
    group.bench_function("evaluate_len_0050", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L50, &emodel))
    });
    group.bench_function("evaluate_len_0100", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L100, &emodel))
    });
    group.bench_function("evaluate_len_0250", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L250, &emodel))
    });
    group.bench_function("evaluate_len_0500", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L500, &emodel))
    });
    group.bench_function("evaluate_len_0750", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L750, &emodel))
    });
    group.bench_function("evaluate_len_1000", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L1000, &emodel))
    });
    group.bench_function("evaluate_len_2500", |b| {
        b.iter(|| evaluate_all_from_file(INPUT_L2500, &emodel))
    });
    group.finish();
}

fn evaluate_all_from_file(path: &str, emodel: &ViennaRNA) {
    let file = File::open(path).expect("Cannot open input file");
    let reader = BufReader::new(file);

    let mut sum = 0i32;
    let mut lines = reader.lines();
    while let Some(Ok(header)) = lines.next() {
        if !header.starts_with('>') {
            panic!("Malformed benchmarking input.");
        }

        let sequence = lines.next().unwrap().unwrap();
        let structure = lines.next().unwrap().unwrap();

        let sequence = NucleotideVec::try_from(sequence.as_str()).unwrap();
        let pairings = PairTable::try_from(structure.as_str())
            .expect("invalid structure in input");
        sum += emodel.energy_of_structure(&sequence, &pairings);
        black_box(sum);
   }
}

criterion_group!(
    benches,
    bench_model_init,
    bench_evaluation,
);
criterion_main!(benches);

