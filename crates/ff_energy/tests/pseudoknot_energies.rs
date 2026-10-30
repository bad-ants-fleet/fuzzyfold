//! Regression test: the energies of the pseudoknotted benchmark sets
//! (lengths 50 and 100) at the default parameters, as one strand, four
//! strands, circular, and four strands with the third circular, as recorded
//! by `benchmarks/fuzzyfold_energies.py` (`benches/data/*_fuzzyfold.txt`).
//! The energies of circular strands with pseudoknots must not depend on
//! where the circle is cut.

use ff_structure::{MultiPairTable, MultiStruct};
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};

const SETS: [(&str, &str); 2] = [
    (include_str!("../benches/data/benchmark_random_pseudoknots_len50.vrna"),
     include_str!("../benches/data/benchmark_random_pseudoknots_len50_fuzzyfold.txt")),
    (include_str!("../benches/data/benchmark_random_pseudoknots_len100.vrna"),
     include_str!("../benches/data/benchmark_random_pseudoknots_len100_fuzzyfold.txt")),
];

/// The string cut into four at the quarters, with the given strand ends.
fn quarters(s: &str, ends: [char; 3]) -> String {
    let n = s.len();
    let cuts = [0, n / 4, n / 2, 3 * n / 4, n];
    let mut out = String::with_capacity(n + 3);
    for k in 0..4 {
        out.push_str(&s[cuts[k]..cuts[k + 1]]);
        if k < 3 {
            out.push(ends[k]);
        }
    }
    out
}

/// One strand, four strands, circular, four strands with the third circular.
fn variants(s: &str) -> [String; 4] {
    [s.to_string(), quarters(s, ['+', '+', '+']), format!("{s}^"), quarters(s, ['+', '+', '^'])]
}

/// The sequences and structures of the sets, with the recorded energies of
/// their variants.
fn records() -> Vec<(&'static str, &'static str, Vec<i32>)> {
    let mut out = Vec::new();
    for (set, energies) in SETS {
        let records: Vec<&str> = set.lines().collect();
        let energies: Vec<&str> = energies.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(records.len(), 3 * energies.len());
        for (record, line) in records.chunks(3).zip(energies) {
            let expected = line.split_whitespace().map(|e| e.parse().unwrap()).collect();
            out.push((record[1], record[2], expected));
        }
    }
    out
}

#[test]
fn pseudoknot_energies_as_recorded() {
    let model = Fuzzyfold::default();
    let records = records();
    assert_eq!(records.len(), 200);
    for (seq, db, expected) in records {
        let found: Vec<i32> = variants(seq).iter().zip(variants(db)).map(|(s, d)| {
            model.energy_of_structure(&NucleotideVec::try_from(s.as_str()).unwrap(),
                &MultiPairTable::try_from(d.as_str()).unwrap()).unwrap()
        }).collect();
        assert_eq!(found, expected, "one four circular four-circular\n{seq}\n{db}");
    }
}

/// The circle cut r bases later: sequence and pair table.
fn rotated(seq: &str, mpt: &MultiPairTable, r: usize) -> (String, MultiPairTable) {
    let m = mpt.len() - 1;
    let bases: Vec<char> = seq.chars().take(m).collect();
    let rseq: String = (0..m).map(|k| bases[(k + r) % m]).chain(['^']).collect();
    let mut rmpt = MultiPairTable::try_from(format!("{}^", ".".repeat(m)).as_str()).unwrap();
    for k in 0..m {
        if let MultiStruct::Paired(p) = mpt.get(k) {
            rmpt.set((k + m - r) % m, MultiStruct::Paired((p + m - r) % m));
        }
    }
    (rseq, rmpt)
}

#[test]
fn circular_pseudoknots_under_rotation() {
    let model = Fuzzyfold::default();
    let energy = |seq: &str, mpt: &MultiPairTable| model.energy_of_structure(
        &NucleotideVec::try_from(seq).unwrap(), mpt).map_err(|e| format!("{e:?}"));
    for (seq, db, _) in records() {
        let (seq, db) = (format!("{seq}^"), format!("{db}^"));
        let mpt = MultiPairTable::try_from(db.as_str()).unwrap();
        let e0 = energy(&seq, &mpt);
        assert!(e0.is_ok(), "{seq}\n{db}");
        for r in 1..mpt.len() - 1 {
            let (rseq, rmpt) = rotated(&seq, &mpt, r);
            assert_eq!(energy(&rseq, &rmpt), e0, "{seq}\n{db} rotated by {r}");
        }
    }
}
