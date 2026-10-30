//! The structures of many trajectories at common output times, and the
//! .drf file format: one line per structure and output time, with an id,
//! the time, the occupancy (the fraction of trajectories), the structure
//! and its energy in kcal/mol; with A-to-I edits also the sequence (as far
//! as transcribed).
//!
//! A structure keeps its id when the strand grows (the id is that of the
//! structure padded to full length), not when it is edited.

use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use rustc_hash::FxHashMap;

use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;

/// A structure with its sequence (full length, with edits) and energy.
pub type State = (NucleotideVec, DotBracketVec, i32);

/// The number of trajectories in each state at each output time.
///
/// ```
/// use std::sync::Arc;
/// use ff_structure::MultiPairTable;
/// use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
/// use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, OutputTimes, SSA, extension_schedule};
/// use ff_kinetics::ensemble::Ensemble;
///
/// let model = Arc::new(Fuzzyfold::default());
/// let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
/// let start = MultiPairTable::try_from("...........").unwrap();
/// let rates = Arrhenius::new(model.temperature(), 1e5);
/// let moves = LoopNeighbors::try_from((seq, &start, model, rates.move_types())).unwrap();
///
/// // 50 times evenly up to 0.1 ms, then 50 evenly in log space up to 1 s.
/// let output = OutputTimes::default();
/// let schedule = extension_schedule(None, 0, output.end()).unwrap();
///
/// let mut ensemble = Ensemble::new(output.times());
/// for _ in 0..10 {
///     let mut ssa = SSA::from((moves.clone(), rates));
///     ensemble.add_trajectory(ssa.sample(&mut rand::rng(), &schedule, &output,
///         |m| (m.current_sequence().clone(), m.current_structure(), m.current_energy())));
/// }
/// // ensemble.write_drf("hairpin.drf", false) writes the file.
/// assert_eq!(ensemble.num_trajectories(), 10);
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble {
    times: Vec<f64>,
    /// `counts[t]`: the number of trajectories in each state at `times[t]`.
    counts: Vec<FxHashMap<State, usize>>,
    num_trajectories: usize,
}

impl Ensemble {
    /// No trajectories at the given output times.
    pub fn new(times: &[f64]) -> Self {
        Self { times: times.to_vec(), counts: vec![FxHashMap::default(); times.len()], num_trajectories: 0 }
    }

    /// Add a trajectory: its state at each output time.
    pub fn add_trajectory(&mut self, states: Vec<State>) {
        assert_eq!(states.len(), self.times.len(), "one state per output time");
        for (counts, state) in self.counts.iter_mut().zip(states) {
            *counts.entry(state).or_insert(0) += 1;
        }
        self.num_trajectories += 1;
    }

    pub fn times(&self) -> &[f64] {
        &self.times
    }

    /// The number of trajectories in each state, per output time.
    pub fn counts(&self) -> &[FxHashMap<State, usize>] {
        &self.counts
    }

    pub fn num_trajectories(&self) -> usize {
        self.num_trajectories
    }

    /// Write a .drf file; with `sequence_column`, the sequence of each
    /// structure is the last column. Per time, the states are sorted by
    /// energy, then structure and sequence; ids are given in that order.
    pub fn write_drf(&self, path: impl AsRef<Path>, sequence_column: bool) -> io::Result<()> {
        let seqlen = self.counts.iter().flat_map(|c| c.keys()).map(|(seq, _, _)| seq.len()).max().unwrap_or(0);
        let id_key = |(seq, s, _): &State| {
            format!("{s}{} {seq}", ".".repeat(seqlen - s.len()))
        };
        let sorted: Vec<Vec<(&State, usize)>> = self.counts.iter().map(|counts| {
            let mut entries: Vec<(&State, usize)> = counts.iter().map(|(s, &c)| (s, c)).collect();
            entries.sort_by_cached_key(|((seq, s, en), _)| (*en, s.to_string(), seq.to_string()));
            entries
        }).collect();
        let mut ids: FxHashMap<String, usize> = FxHashMap::default();
        for (state, _) in sorted.iter().flatten() {
            let next = ids.len();
            ids.entry(id_key(state)).or_insert(next);
        }

        let mut out = BufWriter::new(fs::File::create(path)?);
        write!(out, "id time occupancy structure energy")?;
        writeln!(out, "{}", if sequence_column { " sequence" } else { "" })?;
        for (t, entries) in self.times.iter().zip(&sorted) {
            for (state @ (seq, s, en), count) in entries {
                let occupancy = *count as f64 / self.num_trajectories as f64;
                write!(out, "{:5} {t} {occupancy:5} {s} {:6.2}", ids[&id_key(state)], *en as f64 / 100.0)?;
                if sequence_column {
                    write!(out, " {}", &seq.to_string()[..s.len()])?;
                }
                writeln!(out)?;
            }
        }
        out.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state(seq: &str, s: &str, en: i32) -> State {
        (NucleotideVec::try_from(seq).unwrap(), DotBracketVec::try_from(s).unwrap(), en)
    }

    #[test]
    fn drf_ids_times_and_order() {
        let mut ensemble = Ensemble::new(&[0.0, 2e-6, 1.0]);
        let seq = "GGGGAAACCCC";
        ensemble.add_trajectory(vec![
            state(seq, "....", 0), state(seq, ".....", 0), state(seq, "((((...))))", -443)]);
        ensemble.add_trajectory(vec![
            state(seq, "....", 0), state(seq, ".....", 0), state(seq, "...........", 0)]);
        let path = std::env::temp_dir().join(format!("ff_drf_test_{}.drf", std::process::id()));
        ensemble.write_drf(&path, false).unwrap();
        let text = fs::read_to_string(&path).unwrap();
        fs::remove_file(&path).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "id time occupancy structure energy");
        // The open structure keeps its id while the strand grows; the
        // times are distinct; the lowest energy comes first.
        assert_eq!(lines[1], "    0 0     1 ....   0.00");
        assert_eq!(lines[2], "    0 0.000002     1 .....   0.00");
        assert_eq!(lines[3], "    1 1   0.5 ((((...))))  -4.43");
        assert_eq!(lines[4], "    0 1   0.5 ...........   0.00");
    }
}
