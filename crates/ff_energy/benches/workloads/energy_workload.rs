//! A fixed workload for instruction counts: evaluate a benchmark set of
//! length 1000 in one variant, 20 times, and print the sum of the energies
//! (identical results are part of every comparison).
//!
//! cargo build --release -p ff_energy --example energy_workload
//! valgrind --tool=callgrind target/release/examples/energy_workload structures one
//!
//! Sets: structures, pseudoknots. Variants: one, four, circular,
//! four-circular (four strands, the third circular).

use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/benches/data");

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

fn variant(s: &str, name: &str) -> String {
    match name {
        "one" => s.to_string(),
        "four" => quarters(s, ['+', '+', '+']),
        "circular" => format!("{s}^"),
        "four-circular" => quarters(s, ['+', '+', '^']),
        _ => panic!("Unknown variant: {name}"),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let [_, set, name] = &args[..] else {
        panic!("Usage: energy_workload <structures|pseudoknots> <one|four|circular|four-circular>");
    };
    let input = std::fs::read_to_string(format!("{DATA}/benchmark_random_{set}_len1000.vrna")).unwrap();
    let lines: Vec<&str> = input.lines().collect();
    let inputs: Vec<(NucleotideVec, MultiPairTable)> = lines.chunks(3).map(|c| (
        NucleotideVec::try_from(variant(c[1], name).as_str()).unwrap(),
        MultiPairTable::try_from(variant(c[2], name).as_str()).unwrap(),
    )).collect();
    let model = Fuzzyfold::default();
    let mut sum: i64 = 0;
    for _ in 0..20 {
        for (seq, mpt) in &inputs {
            sum += model.energy_of_structure(seq, mpt).unwrap() as i64;
        }
    }
    println!("{set} {name} {sum}");
}
