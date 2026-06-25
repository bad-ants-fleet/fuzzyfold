use std::fs::File;
use std::sync::Arc;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Cursor;
use std::hint::black_box;
use criterion::BatchSize;
use criterion::criterion_group;
use criterion::criterion_main;
use criterion::Criterion;
use rand::SeedableRng;
use rand::rngs::StdRng;

use ff_structure::PairTable;
use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA;
use ff_kinetics::Walker;
use ff_kinetics::Arrhenius;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy;
use ff_kinetics::SSA;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::MacrostateRegistryPL;

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

const MACROSTATES_L50: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len50.ms");
const MACROSTATES_L100: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len100.ms");
const MACROSTATES_L250: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len250.ms");
const MACROSTATES_L500: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len500.ms");
const MACROSTATES_L750: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len750.ms");
const MACROSTATES_L1000: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len1000.ms");
const MACROSTATES_L2500: &str = concat!(env!("CARGO_MANIFEST_DIR"),
    "/benches/data/benchmark_macrostates_len2500.ms");

struct BenchCase {
    name: &'static str,
    inputs_path: &'static str,
    macrostates_path: &'static str,
}

const CASES: &[BenchCase] = &[
    BenchCase { name: "len_0050", inputs_path: INPUT_L50,  macrostates_path: MACROSTATES_L50  },
    BenchCase { name: "len_0100", inputs_path: INPUT_L100, macrostates_path: MACROSTATES_L100 },
    BenchCase { name: "len_0250", inputs_path: INPUT_L250, macrostates_path: MACROSTATES_L250 },
    BenchCase { name: "len_0500", inputs_path: INPUT_L500, macrostates_path: MACROSTATES_L500 },
    BenchCase { name: "len_0750", inputs_path: INPUT_L750, macrostates_path: MACROSTATES_L750 },
    BenchCase { name: "len_1000", inputs_path: INPUT_L1000, macrostates_path: MACROSTATES_L1000 },
    BenchCase { name: "len_2500", inputs_path: INPUT_L2500, macrostates_path: MACROSTATES_L2500 },
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

/// Split a multi-macrostate file on '>' lines, returning each block as a String.
fn split_macrostate_blocks(path: &str) -> Vec<String> {
    let file = File::open(path).expect("Cannot open macrostate file");
    let reader = BufReader::new(file);
    let mut blocks: Vec<String> = Vec::new();
    let mut current = String::new();

    for line in reader.lines() {
        let line = line.unwrap();
        if line.trim_start().starts_with('>') && !current.is_empty() {
            blocks.push(current.trim_end().to_string());
            current = String::new();
        }
        current.push_str(&line);
        current.push('\n');
    }
    if !current.trim().is_empty() {
        blocks.push(current.trim_end().to_string());
    }
    blocks
}

fn load_dbv_registries(
    macrostates_path: &str,
    inputs: &[(NucleotideVec, PairTable)],
    emodel: &Arc<ViennaRNA>,
) -> Vec<MacrostateRegistry<ViennaRNA>> {
    let blocks = split_macrostate_blocks(macrostates_path);
    assert_eq!(blocks.len(), inputs.len(),
        "Macrostate block count must match input trajectory count");
    blocks.iter().zip(inputs.iter()).map(|(block, (seq, _))| {
        let mut registry = MacrostateRegistry::from((seq.clone(), emodel.clone()));
        registry.insert_from_reader(Cursor::new(block.as_bytes()), macrostates_path)
            .expect("Failed to parse macrostate block");
        registry
    }).collect()
}

fn load_pl_registries(
    macrostates_path: &str,
    inputs: &[(NucleotideVec, PairTable)],
    emodel: &Arc<ViennaRNA>,
) -> Vec<MacrostateRegistryPL<ViennaRNA>> {
    let blocks = split_macrostate_blocks(macrostates_path);
    assert_eq!(blocks.len(), inputs.len(),
        "Macrostate block count must match input trajectory count");
    blocks.iter().zip(inputs.iter()).map(|(block, (seq, _))| {
        let mut registry = MacrostateRegistryPL::from((seq.clone(), emodel.clone()));
        registry.insert_from_reader(Cursor::new(block.as_bytes()), macrostates_path)
            .expect("Failed to parse macrostate block");
        registry
    }).collect()
}

// ---------------------------------------------------------------------------
// DBV registry
// ---------------------------------------------------------------------------

fn bench_classify_dbv(c: &mut Criterion) {
    let emodel = Arc::new(ViennaRNA::default());
    let rmodel = Arrhenius::new(emodel.temperature(), 1.0, None, None);
    let mut group = c.benchmark_group("Classify macrostates (DBV) during SSA");
    group.measurement_time(std::time::Duration::from_secs(50));

    for case in CASES {
        let inputs = load_raw_inputs(case.inputs_path);
        let registries = load_dbv_registries(case.macrostates_path, &inputs, &emodel);
        let mut rng = StdRng::seed_from_u64(42);

        group.bench_function(format!("dbv_{}", case.name), |b| {
            b.iter_batched(
                || inputs.iter().zip(registries.iter()).collect::<Vec<_>>(),
                |pairs| {
                    for ((seq, pt), registry) in pairs {
                        let moves = LoopNeighbors::try_from((
                            seq.clone(), pt, emodel.clone(), shift_policy::NoShift,
                        ))
                        .expect("Failed to build LoopNeighbors");
                        let mut simulator = SSA::from((moves, rmodel));
                        let mut classifications = 0usize;

                        simulator.simulate(
                            &mut rng,
                            black_box(10.0),
                            |_t, _tinc, _rsum, walker| {
                                let s = walker.current_structure();
                                classifications += registry.classify(&s);
                                true
                            },
                        );
                        black_box(classifications);
                    }
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// PL registry
// ---------------------------------------------------------------------------

fn bench_classify_pl(c: &mut Criterion) {
    let emodel = Arc::new(ViennaRNA::default());
    let rmodel = Arrhenius::new(emodel.temperature(), 1.0, None, None);
    let mut group = c.benchmark_group("Classify macrostates (PL) during SSA");
    group.measurement_time(std::time::Duration::from_secs(50));

    for case in CASES {
        let inputs = load_raw_inputs(case.inputs_path);
        let registries = load_pl_registries(case.macrostates_path, &inputs, &emodel);
        let mut rng = StdRng::seed_from_u64(42);

        group.bench_function(format!("pl_{}", case.name), |b| {
            b.iter_batched(
                || inputs.iter().zip(registries.iter()).collect::<Vec<_>>(),
                |pairs| {
                    for ((seq, pt), registry) in pairs {
                        let moves = LoopNeighbors::try_from((
                            seq.clone(), pt, emodel.clone(), shift_policy::NoShift,
                        ))
                        .expect("Failed to build LoopNeighbors");
                        let mut simulator = SSA::from((moves, rmodel));
                        let mut classifications = 0usize;

                        simulator.simulate(
                            &mut rng,
                            black_box(10.0),
                            |_t, _tinc, _rsum, walker| {
                                let s = walker.current_structure();
                                classifications += registry.classify(&s);
                                true
                            },
                        );
                        black_box(classifications);
                    }
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// Isolated classify() — raw lookup cost without simulation noise.
// ---------------------------------------------------------------------------

fn bench_classify_isolated(c: &mut Criterion) {
    let emodel = Arc::new(ViennaRNA::default());
    let mut group = c.benchmark_group("Isolated classify() calls");

    let inputs = load_raw_inputs(INPUT_L50);
    let structures: Vec<DotBracketVec> = inputs.iter()
        .map(|(_, pt)| DotBracketVec::try_from(pt).unwrap())
        .collect();

    let dbv_registries = load_dbv_registries(MACROSTATES_L50, &inputs, &emodel);
    let pl_registries  = load_pl_registries(MACROSTATES_L50, &inputs, &emodel);

    group.bench_function("isolated_dbv", |b| {
        b.iter(|| {
            let mut sum = 0usize;
            for (s, registry) in structures.iter().zip(dbv_registries.iter()) {
                sum += registry.classify(black_box(s));
            }
            black_box(sum)
        })
    });

    group.bench_function("isolated_pl", |b| {
        b.iter(|| {
            let mut sum = 0usize;
            for (s, registry) in structures.iter().zip(pl_registries.iter()) {
                sum += registry.classify(black_box(s));
            }
            black_box(sum)
        })
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_classify_dbv,
    bench_classify_pl,
    bench_classify_isolated,
);
criterion_main!(benches);
