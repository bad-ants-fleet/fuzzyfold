//! Seeded trajectories, with and without shift moves, on linear and
//! circular strands, compared with recorded digests (`data/trajectories.txt`).
//! At every step, the flux is the sum of the rates of the moves of a move set
//! built from scratch, and the energy is that of the structure.
//!
//! To record the digests again: `FF_RECORD_TRAJECTORIES=1 cargo test -p
//! ff_kinetics --test trajectories`.

use std::sync::Arc;
use rand::SeedableRng;
use rand::rngs::StdRng;
use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
use ff_kinetics::{shift_policy, shift_policy::ShiftPolicy, Arrhenius, LoopNeighbors, RateModel, Walker, SSA};

const RUNS: [(&str, &str); 4] = [
    ("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC", ".((((....)))).((((........))))..............."),
    ("GGGAGCUAAAGCUCCCAUCGGCAUUAGCCGAUGGGAAAGCAUCCAGCUUAGG", "...................................................."),
    ("GCGCAAAAGCGCAAAAAGCGCUUUUGCGC^", ".............................^"),
    ("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC^", ".((((....)))).((((........))))...............^"),
];

const STEPS: usize = 1000;

/// FNV-1a, to keep the digest of a trajectory short and exact.
fn digest(lines: &[String]) -> u64 {
    lines.iter().flat_map(|l| l.bytes().chain(*b"\n"))
        .fold(0xcbf29ce484222325, |h, b| (h ^ b as u64).wrapping_mul(0x100000001b3))
}

/// The trajectory of a run, one line per step: time, waiting time, flux,
/// energy and structure, with floats as their exact shortest form.
fn trajectory<P: ShiftPolicy + Copy>(seq: &str, db: &str, policy: P, shifts: bool) -> Vec<String> {
    let model = Arc::new(Fuzzyfold::default());
    let shift = shifts.then_some(1.0);
    let rates = Arrhenius::new(model.temperature(), 1.0, shift, shift);
    let sequence = NucleotideVec::try_from(seq).unwrap();
    let start = MultiPairTable::try_from(db).unwrap();
    let moves = LoopNeighbors::try_from((sequence.clone(), &start, model.clone(), policy)).unwrap();
    let mut ssa = SSA::from((moves, rates));
    let mut lines = Vec::new();
    ssa.simulate(&mut StdRng::seed_from_u64(2026), 1e9, |t, tinc, flux, w| {
        let structure = w.current_structure();
        let pairs = MultiPairTable::try_from(&structure).unwrap();
        assert_eq!(w.current_energy(), model.energy_of_structure(&sequence, &pairs).unwrap(), "{structure}");
        let fresh = LoopNeighbors::try_from((sequence.clone(), &pairs, model.clone(), policy)).unwrap();
        let total: f64 = fresh.propose_moves().map(|(m, delta)| rates.rate(&m, delta)).sum();
        assert!((flux - total).abs() <= 1e-9 * total, "{structure}: flux {flux}, rates {total}");
        lines.push(format!("{t:?} {tinc:?} {flux:?} {} {structure}", w.current_energy()));
        lines.len() < STEPS
    });
    lines
}

#[test]
fn trajectories_as_recorded() {
    let mut records = Vec::new();
    for (seq, db) in RUNS {
        for (name, lines) in [
            ("no-shift", trajectory(seq, db, shift_policy::NoShift, false)),
            ("shift", trajectory(seq, db, shift_policy::ThreeAndFour, true)),
        ] {
            records.push(format!("{seq} {name} {} {:016x} {}", lines.len(), digest(&lines), lines.last().unwrap()));
        }
    }
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/trajectories.txt");
    if std::env::var_os("FF_RECORD_TRAJECTORIES").is_some() {
        std::fs::write(path, records.join("\n") + "\n").unwrap();
        return;
    }
    let recorded = std::fs::read_to_string(path).unwrap();
    for (now, then) in records.iter().zip(recorded.lines()) {
        assert_eq!(now, then);
    }
    assert_eq!(records.len(), recorded.lines().count());
}
