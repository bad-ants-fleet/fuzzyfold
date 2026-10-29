use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::hint::black_box;
use criterion::BatchSize;
use criterion::Criterion;
use criterion::criterion_main;
use criterion::criterion_group;

use ff_structure::MultiPairTable;
use ff_energy::EnergyModel;
use ff_energy::Fuzzyfold;
use ff_energy::NucleotideVec;

const LENGTHS: [usize; 3] = [100, 1000, 2500];

/// The benchmark name of a sequence length, e.g. "len_0100".
fn name(len: usize) -> String {
    format!("len_{len:04}")
}

/// The sequence and structure lines of benchmark_random_{set}_len{len}.vrna.
fn load_lines(set: &str, len: usize) -> Vec<(String, String)> {
    let path = format!("{}/benches/data/benchmark_random_{set}_len{len}.vrna",
        env!("CARGO_MANIFEST_DIR"));
    let file = File::open(path).expect("Cannot open input file");
    let lines: Vec<String> = BufReader::new(file).lines().map(|l| l.unwrap()).collect();
    lines.chunks(3).map(|c| {
        assert!(c[0].starts_with('>'), "Malformed benchmarking input.");
        (c[1].clone(), c[2].clone())
    }).collect()
}

/// The string cut into four at the quarters, with the given strand ends.
fn quarters(s: &str, ends: [char; 3]) -> String {
    let n = s.len();
    let cuts = [n / 4, n / 2, 3 * n / 4];
    let mut out = String::with_capacity(n + 3);
    let mut from = 0;
    for (cut, end) in cuts.into_iter().zip(ends) {
        out.push_str(&s[from..cut]);
        out.push(end);
        from = cut;
    }
    out.push_str(&s[from..]);
    out
}

/// The variants of the benchmark: one strand, four strands, circular, and
/// four strands with the third circular.
const VARIANTS: [&str; 4] = ["one", "four", "circular", "four-circular"];

fn variant(s: &str, name: &str) -> String {
    match name {
        "one" => s.to_string(),
        "four" => quarters(s, ['+', '+', '+']),
        "circular" => format!("{s}^"),
        "four-circular" => quarters(s, ['+', '+', '^']),
        _ => unreachable!(),
    }
}

fn load_inputs(set: &str, len: usize, name: &str) -> Vec<(NucleotideVec, MultiPairTable)> {
    load_lines(set, len).iter().map(|(seq, db)| (
        NucleotideVec::try_from(variant(seq, name).as_str()).unwrap(),
        MultiPairTable::try_from(variant(db, name).as_str()).unwrap(),
    )).collect()
}

fn bench_model_init(c: &mut Criterion) {
    c.bench_function("Fuzzyfold::default()", |b| {
        b.iter(|| {
            let model = Fuzzyfold::default();
            black_box(model);
        })
    });
}

/// The evaluation of a set in each of its variants: the nested structures
/// in one and four strands, the pseudoknotted ones in all variants.
fn bench_set(c: &mut Criterion, set: &str, variants: &[&str]) {
    let emodel = Fuzzyfold::default();
    let mut group = c.benchmark_group(format!("Evaluation of {set}"));
    for len in LENGTHS {
        for &name in variants {
            let inputs = load_inputs(set, len, name);
            group.bench_function(format!("{name}_{}", self::name(len)), |b| {
                b.iter_batched(
                    || &inputs,
                    |inputs| {
                        let mut sum = 0i32;
                        for (seq, pt) in inputs {
                            sum += emodel.energy_of_structure(black_box(seq), black_box(pt)).unwrap();
                        }
                        black_box(sum);
                    },
                    BatchSize::LargeInput,
                )
            });
        }
    }
    group.finish();
}

fn bench_structures(c: &mut Criterion) {
    bench_set(c, "structures", &VARIANTS[..2]);
}

fn bench_pseudoknots(c: &mut Criterion) {
    bench_set(c, "pseudoknots", &VARIANTS);
}

criterion_group!(
    benches,
    bench_model_init,
    bench_structures,
    bench_pseudoknots,
);
criterion_main!(benches);
