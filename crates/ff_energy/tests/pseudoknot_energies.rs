//! Regression test: the energies of random structures with pseudoknots, on
//! one strand and on several, linear and circular, at the default
//! parameters (`data/pseudoknot_energies.txt`). The energies of circular
//! strands with pseudoknots must not depend on where the circle is cut;
//! those recorded as errors are not checked.

use ff_structure::{MultiPairTable, MultiStruct};
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};

fn records() -> Vec<(&'static str, &'static str, &'static str)> {
    include_str!("data/pseudoknot_energies.txt").lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| match l.split_whitespace().collect::<Vec<_>>()[..] {
            [seq, db, e] => (seq, db, e),
            _ => panic!("{l}"),
        })
        .collect()
}

#[test]
fn pseudoknot_energies_as_recorded() {
    let model = Fuzzyfold::default();
    let records = records();
    assert_eq!(records.len(), 1500);
    for (seq, db, expected) in records {
        if db.ends_with('^') && expected == "error" {
            continue;
        }
        let mpt = MultiPairTable::try_from(db).unwrap();
        let energy = model.energy_of_structure(&NucleotideVec::try_from(seq).unwrap(), &mpt);
        match expected {
            "error" => assert!(energy.is_err(), "{seq}\n{db} {energy:?}"),
            e => assert_eq!(energy.ok().map(|en| format!("{:.2}", en as f64 / 100.)).as_deref(), Some(e), "{seq}\n{db}"),
        }
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
    let mut tested = 0;
    for (seq, db, _) in records().into_iter().filter(|r| r.1.ends_with('^')) {
        let mpt = MultiPairTable::try_from(db).unwrap();
        let e0 = energy(seq, &mpt);
        for r in 1..mpt.len() - 1 {
            let (rseq, rmpt) = rotated(seq, &mpt, r);
            assert_eq!(energy(&rseq, &rmpt), e0, "{seq}\n{db} rotated by {r}");
        }
        tested += e0.is_ok() as usize;
    }
    assert!(tested > 50, "{tested} circles evaluated");
}
