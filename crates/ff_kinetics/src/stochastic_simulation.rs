use ff_structure::DotBracketVec;
use rand::{Rng, RngExt};

use crate::MoveSet;
use crate::Elongation;
use crate::Neighbors;
use crate::RateModel;
use crate::rate_tree::RateTree;

/// The stochastic simulation algorithm (Gillespie) on a move set: the rates
/// of all moves of the current structure in a tree, for drawing the next
/// move in logarithmic time.
pub struct SSA<M: MoveSet, K: RateModel> {
    moves: M,
    ratemodel: K,
    rate_tree: RateTree,
}

impl<M: MoveSet, K: RateModel> From<(M, K)> for SSA<M, K> {
    fn from((moves, ratemodel): (M, K)) -> Self {
        let rate_tree = moves.propose_moves()
            .map(|(mv, delta)| (mv, ratemodel.rate(&mv, delta)))
            .filter(|&(_, k)| k > 0.0)
            .collect();
        Self { moves, ratemodel, rate_tree }
    }
}

impl<M: MoveSet, K: RateModel> SSA<M, K> {
    pub fn move_set(&self) -> &M {
        &self.moves
    }

    pub fn current_structure(&self) -> DotBracketVec {
        self.moves.current_structure()
    }

    pub fn current_energy(&self) -> i32 {
        self.moves.current_energy()
    }

    /// A waiting time in the current structure, drawn from the exponential
    /// distribution of the total rate (the flux), and the flux. None if no
    /// move has a rate.
    pub fn waiting_time<R: Rng + ?Sized>(&self, rng: &mut R) -> Option<(f64, f64)> {
        let flux = self.rate_tree.total_rate();
        (flux > 0.0).then(|| (-rng.random::<f64>().ln() / flux, flux))
    }

    /// Apply a move, drawn with a probability proportional to its rate.
    /// Panics if no move has a rate.
    pub fn apply_random_move<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let threshold = rng.random::<f64>() * self.rate_tree.total_rate();
        let mv = self.rate_tree.select_by_threshold(threshold).expect("Must select a move!");
        let (old, new) = self.moves.apply_move(&mv);
        self.update_rate_tree(old, new);
    }

    /// Simulate for the time t_max. For each structure, the callback gets
    /// the time, the waiting time, the flux and the move set; the
    /// simulation stops when it returns false. A move whose time falls
    /// after t_max is not applied: the simulation ends in the structure of
    /// the last callback. Without moves, the callback is called once with
    /// the rest of the time. Returns the last return value of the callback.
    ///
    /// ```
    /// use std::sync::Arc;
    /// use ff_structure::{DotBracketVec, MultiPairTable};
    /// use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
    /// use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, SSA};
    ///
    /// let model = Arc::new(Fuzzyfold::default());
    /// let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
    /// let start = MultiPairTable::try_from("...........").unwrap();
    /// let rates = Arrhenius::new(model.temperature(), 1e5);
    /// let moves = LoopNeighbors::try_from((seq, &start, model, rates.move_types())).unwrap();
    /// let mut ssa = SSA::from((moves, rates));
    ///
    /// // The first passage time to the hairpin, within at most 1 s.
    /// let hairpin = DotBracketVec::try_from("((((...))))").unwrap();
    /// let mut first_passage = None;
    /// ssa.simulate(&mut rand::rng(), 1.0, |t, _waiting_time, _flux, m| {
    ///     if m.current_structure() == hairpin {
    ///         first_passage = Some(t);
    ///     }
    ///     first_passage.is_none()
    /// });
    /// println!("{first_passage:?}");
    /// ```
    pub fn simulate<R, F>(&mut self, rng: &mut R, t_max: f64, mut callback: F) -> bool
    where
        R: Rng + ?Sized,
        F: FnMut(f64, f64, f64, &M) -> bool,
    {
        let mut t = 0.;
        while t < t_max {
            let Some((tinc, flux)) = self.waiting_time(rng) else {
                return callback(t, t_max - t, 0.0, &self.moves);
            };
            if !callback(t, tinc, flux, &self.moves) {
                return false;
            }
            t += tinc;
            if t < t_max {
                self.apply_random_move(rng);
            }
        }
        true
    }

    /// Apply the changes of a move (see `MoveSet::apply_move`): a new move
    /// takes the place of an existing one with the same key, else of an old
    /// one, else it is inserted; old moves without a successor are removed.
    /// Moves with rate 0 are not in the tree: their rate constant is 0, so
    /// a move never changes between rate 0 and a positive rate.
    fn update_rate_tree(&mut self, old: Neighbors, new: Neighbors) {
        let mut del = old.iter();
        let mut add = new.iter();
        let mut cur_del = del.next();
        let mut cur_add = add.next();
        while cur_del.is_some() || cur_add.is_some() {
            match (cur_del, cur_add) {
                (Some((omv, _)), Some((nmv, delta))) => {
                    let k = self.ratemodel.rate(nmv, *delta);
                    if k == 0.0 || self.rate_tree.update_rate(nmv, k) {
                        cur_add = add.next();
                    } else if self.rate_tree.replace(omv, nmv, k) {
                        // Only true if the old move exists: moves with
                        // rate 0 are not in the tree.
                        cur_add = add.next();
                        cur_del = del.next();
                    } else {
                        cur_del = del.next();
                    }
                }
                (Some((omv, _)), None) => {
                    self.rate_tree.remove(omv);
                    cur_del = del.next();
                }
                (None, Some((nmv, delta))) => {
                    let k = self.ratemodel.rate(nmv, *delta);
                    if k > 0. && !self.rate_tree.update_rate(nmv, k) {
                        self.rate_tree.insert(*nmv, k);
                    }
                    cur_add = add.next();
                }
                (None, None) => unreachable!(),
            }
        }
    }
}

impl<M: Elongation, K: RateModel> SSA<M, K> {
    /// Extend the structure by one unpaired base at its 3' end.
    pub fn apply_extension(&mut self) {
        let (old, new) = self.moves.apply_extension();
        self.update_rate_tree(old, new);
    }

    /// Simulate transcription: for `times[k]` at each length, from the length
    /// of the current structure to that of the sequence, then extend. The
    /// callback is as for `simulate`, with the time since the start, and a
    /// waiting time that ends at most slightly after the next extension.
    ///
    /// ```
    /// use std::sync::Arc;
    /// use ff_structure::MultiPairTable;
    /// use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
    /// use ff_kinetics::{Arrhenius, LoopNeighbors, SSA, extension_schedule};
    ///
    /// let model = Arc::new(Fuzzyfold::default());
    /// let seq = NucleotideVec::try_from_rna("GGGGAAACCCCAAAGGGGAAACCCC").unwrap();
    /// // Transcription starts with the first base.
    /// let start = MultiPairTable::try_from(".").unwrap();
    /// let rates = Arrhenius::new(model.temperature(), 1e5);
    /// let moves = LoopNeighbors::try_from((seq.clone(), &start, model, rates.move_types())).unwrap();
    /// let mut ssa = SSA::from((moves, rates));
    ///
    /// // One base every 20 ms, then 1 s at full length.
    /// let num_ext = seq.len() - 1;
    /// let schedule = extension_schedule(Some(0.02), num_ext, 0.02 * num_ext as f64 + 1.0).unwrap();
    /// ssa.co_simulate(&mut rand::rng(), &schedule, |_t, _waiting_time, _flux, _m| true);
    /// println!("{} {}", ssa.current_structure(), ssa.current_energy());
    /// ```
    pub fn co_simulate<R, F>(&mut self, rng: &mut R, times: &[f64], mut callback: F)
    where
        R: Rng + ?Sized,
        F: FnMut(f64, f64, f64, &M) -> bool,
    {
        //NOTE: This panic is temporary, should return an error eventually.
        if self.moves.sequence_length() != self.moves.current_structure_length() - 1 + times.len() {
            panic!("mismatch between simulation times and sequence length");
        }

        let mut gtime = 0.0;
        for (idx, &time) in times.iter().enumerate() {
            // The waiting time reported to the callback ends at most eps
            // after the extension, so that a caller that records structures
            // at output times needs no check of its own against the
            // extension times.
            let eps = time / 1000.0;
            let co_callback = |t: f64, tinc: f64, flux: f64, m: &M| {
                callback(t + gtime, tinc.min(time - t + eps), flux, m)
            };
            let cb = self.simulate(rng, time, co_callback);

            // Skip extension after the last time point
            if !cb || idx + 1 == times.len() {
                break;
            }
            self.apply_extension();
            gtime += time;
        }
    }

    /// Simulate transcription as `co_simulate` (see `extension_schedule`)
    /// and record the move set at each output time. The schedule must end
    /// at the last output time.
    ///
    /// ```
    /// use std::sync::Arc;
    /// use ff_structure::MultiPairTable;
    /// use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
    /// use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, OutputTimes, SSA, extension_schedule};
    ///
    /// let model = Arc::new(Fuzzyfold::default());
    /// let seq = NucleotideVec::try_from_rna("GGGGAAACCCCAAAGGGGAAACCCC").unwrap();
    /// let start = MultiPairTable::try_from(".").unwrap();
    /// let rates = Arrhenius::new(model.temperature(), 1e5);
    /// let moves = LoopNeighbors::try_from((seq.clone(), &start, model, rates.move_types())).unwrap();
    /// let mut ssa = SSA::from((moves, rates));
    ///
    /// // One base every 20 ms: the structure just before each extension,
    /// // then at 50 times up to 1 s after the full length.
    /// let num_ext = seq.len() - 1;
    /// let output = OutputTimes::transcription(0.02, num_ext, 1.0).unwrap();
    /// let schedule = extension_schedule(Some(0.02), num_ext, output.end()).unwrap();
    /// let structures = ssa.sample(&mut rand::rng(), &schedule, &output, |m| m.current_structure());
    /// for (t, s) in output.times().iter().zip(&structures) {
    ///     println!("{t:.3} {s}");
    /// }
    /// ```
    pub fn sample<R, T, F>(&mut self, rng: &mut R, schedule: &[f64], output: &OutputTimes,
        mut record: F) -> Vec<T>
    where
        R: Rng + ?Sized,
        F: FnMut(&M) -> T,
    {
        let times = output.times();
        let mut samples = Vec::with_capacity(times.len());
        self.co_simulate(rng, schedule, |t, tinc, _, m| {
            while samples.len() < times.len() && t + tinc >= times[samples.len()] {
                samples.push(record(m));
            }
            true
        });
        // Output times at the very end, after the last reported waiting time.
        while samples.len() < times.len() {
            samples.push(record(&self.moves));
        }
        samples
    }
}

/// The times at which a simulation reports its state, from 0.
///
/// ```
/// use ff_kinetics::OutputTimes;
///
/// // 50 times evenly up to 0.1 ms, then 50 evenly in log space up to 1 s.
/// let default = OutputTimes::default();
/// assert_eq!(default.end(), 1.0);
///
/// // 10 times evenly up to 1 ms, then 90 in log space up to 10 s.
/// let longer = OutputTimes::lin_log(1e-3, 10.0, 10, 90).unwrap();
/// assert_eq!(longer.t_sep(), Some(1e-3));
///
/// // Any increasing list.
/// let three = OutputTimes::new(vec![0.0, 0.5, 1.0]).unwrap();
/// assert_eq!(three.times().len(), 3);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct OutputTimes {
    times: Vec<f64>,
    t_sep: Option<f64>,
}

impl OutputTimes {
    /// Any increasing list of times, from 0 or later.
    pub fn new(times: Vec<f64>) -> Result<Self, String> {
        if times.is_empty() {
            return Err("no output times".to_string());
        }
        if times[0] < 0.0 || times.windows(2).any(|w| w[0] >= w[1]) {
            return Err("output times must increase from 0 or later".to_string());
        }
        Ok(Self { times, t_sep: None })
    }

    /// 0, then t_lin times evenly up to t_sep, then t_log times evenly in
    /// log space up to t_end.
    pub fn lin_log(t_sep: f64, t_end: f64, t_lin: usize, t_log: usize) -> Result<Self, String> {
        if !(0.0 < t_sep && t_sep < t_end) {
            return Err(format!("output times need 0 < t_sep ({t_sep}) < t_end ({t_end})"));
        }
        if t_lin == 0 || t_log == 0 {
            return Err("output times need t_lin > 0 and t_log > 0".to_string());
        }
        let step = t_sep / t_lin as f64;
        let mut times = vec![0.0];
        times.extend((1..=t_lin).map(|i| i as f64 * step));
        let (log_sep, log_end) = (times[t_lin].ln(), t_end.ln());
        times.extend((1..t_log).map(|i| (log_sep + i as f64 / t_log as f64 * (log_end - log_sep)).exp()));
        times.push(t_end);
        Ok(Self { times, t_sep: Some(t_sep) })
    }

    /// Co-transcriptional output times: one per extension (every t_ext) up
    /// to the full length, then 50 in log space up to t_end after it.
    pub fn transcription(t_ext: f64, num_ext: usize, t_end: f64) -> Result<Self, String> {
        if num_ext == 0 {
            return Err("co-transcriptional output times need extensions".to_string());
        }
        let t_sep = t_ext * num_ext as f64;
        Self::lin_log(t_sep, t_sep + t_end, num_ext, 50)
    }

    pub fn times(&self) -> &[f64] {
        &self.times
    }

    /// The end of the linear part, if made by lin_log.
    pub fn t_sep(&self) -> Option<f64> {
        self.t_sep
    }

    /// The last output time, where a simulation ends.
    pub fn end(&self) -> f64 {
        *self.times.last().expect("output times are not empty")
    }
}

/// 50 times evenly up to 1e-4 s, then 50 in log space up to 1 s.
impl Default for OutputTimes {
    fn default() -> Self {
        Self::lin_log(1e-4, 1.0, 50, 50).expect("valid default output times")
    }
}

/// The times of `co_simulate`: t_ext at each length before the full length
/// (num_ext extensions), then the rest up to the end.
pub fn extension_schedule(t_ext: Option<f64>, num_ext: usize, end: f64) -> Result<Vec<f64>, String> {
    match t_ext {
        None if num_ext > 0 => Err("t_ext is needed when the start is shorter than the sequence".to_string()),
        None => Ok(vec![end]),
        Some(_) if num_ext == 0 => Err("t_ext is given, but the start has the length of the sequence".to_string()),
        Some(t_ext) if t_ext <= 0.0 => Err(format!("t_ext ({t_ext}) must be positive")),
        Some(t_ext) => {
            let transcription = t_ext * num_ext as f64;
            if end <= transcription {
                return Err(format!("the simulation ends ({end}) before the full length ({transcription})"));
            }
            let mut schedule = vec![t_ext; num_ext];
            schedule.push(end - transcription);
            Ok(schedule)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rand::rng;
    use std::sync::Arc;
    use std::collections::HashSet;

    use ff_structure::MultiPairTable;
    use ff_energy::Fuzzyfold;
    use ff_energy::EnergyModel;
    use ff_energy::NucleotideVec;
    use crate::Arrhenius;
    use crate::movesets::LoopNeighbors;
    use crate::MoveTypes;
    use crate::LoopTable;
    use crate::movesets::Move;
    use rustc_hash::FxHashMap;

    macro_rules! setup_ssa_input {
        ($wname:ident, $rname:ident, $seq:expr, $db:expr) => {
            let emodel = Fuzzyfold::default();
            let $rname = Arrhenius::new(emodel.temperature(), 1.0);

            let sequence = NucleotideVec::try_from($seq)
                .expect("Invalid sequence?");
            let pairings = MultiPairTable::try_from($db)
                .expect("Invalid structure?");

            let ltab = LoopTable::try_from((sequence, &pairings, Arc::new(emodel)))
                .expect("Invalid sequence/structure combination");
            let $wname = LoopNeighbors::from((ltab, MoveTypes::NONE));
        };
    }

    fn same_moves(a: &[(Move, i32)], b: &[(Move, i32)]) -> bool {
        let sa: HashSet<_> = a.iter().copied().collect();
        let sb: HashSet<_> = b.iter().copied().collect();
        sa == sb
    }

    /// The moves and rates in the tree are those of a move set built from
    /// scratch on the current sequence, after every move and extension,
    /// with shift moves and A-to-I edits, also co-transcriptionally and on
    /// a circle.
    #[test]
    fn rate_tree_after_moves_and_extensions() {
        let emodel = Arc::new(Fuzzyfold::default());
        let rmodel = Arrhenius::new(emodel.temperature(), 1.0)
            .with_three_way(1.0).with_four_way(1.0).with_a_to_i(1.0);
        let types = rmodel.move_types();
        let runs = [
            ("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC", "."),
            ("GGGAGCUAAAGCUCCCAUCGGCAUUAGCCGAUGGGAAAGCAUCC", "((((((....))))))"),
            ("GCGCAAAAGCGCAAAAAGCGCUUUUGCGC^", ".............................^"),
            ("GCGACGCAAAAGCGCCGCAAGCGACGCAAAAGCGCCGC", "(((.(((....))).)))..(((.(((....))).)))"),
        ];
        let mut rng = StdRng::seed_from_u64(3);
        let mut edited = 0;
        for (seq, db) in runs {
            let sequence = NucleotideVec::try_from(seq).unwrap();
            let start = MultiPairTable::try_from(db).unwrap();
            let moves = LoopNeighbors::try_from((sequence.clone(), &start, emodel.clone(), types)).unwrap();
            let mut ssa = SSA::from((moves, rmodel));
            for step in 0..3000 {
                let length = ssa.move_set().current_structure_length();
                let stuck = ssa.rate_tree.is_empty();
                if length < sequence.len() && (stuck || step % 20 == 19) {
                    ssa.apply_extension();
                } else if stuck {
                    break;
                } else {
                    ssa.apply_random_move(&mut rng);
                }
                let db = ssa.current_structure();
                let pairs = MultiPairTable::try_from(&db).unwrap();
                let current = ssa.move_set().loop_table().sequence().clone();
                let fresh = LoopNeighbors::try_from((current, &pairs, emodel.clone(), types)).unwrap();
                let expected: FxHashMap<Move, f64> = fresh.propose_moves()
                    .map(|(mv, delta)| (mv, rmodel.rate(&mv, delta)))
                    .filter(|&(_, k)| k > 0.0)
                    .collect();
                let tree: FxHashMap<Move, f64> = ssa.rate_tree.iter().collect();
                assert_eq!(tree, expected, "{seq} step {step}: {db}");
            }
            assert_eq!(ssa.move_set().current_structure_length(), sequence.len());
            edited += ssa.move_set().loop_table().sequence().iter().filter(|&&b| b == ff_energy::Base::IN).count();
        }
        assert!(edited > 0, "no edit happened");
    }

    #[test]
    fn output_times_and_schedule() {
        let out = OutputTimes::lin_log(1e-3, 10.0, 4, 3).unwrap();
        assert_eq!(out.times().len(), 1 + 4 + 3);
        assert_eq!(out.times()[..5], [0.0, 2.5e-4, 5e-4, 7.5e-4, 1e-3]);
        assert_eq!((out.t_sep(), out.end()), (Some(1e-3), 10.0));
        assert!(out.times().windows(2).all(|w| w[0] < w[1]));
        assert_eq!(OutputTimes::default().end(), 1.0);
        assert!(OutputTimes::new(vec![0.0, 1.0, 1.0]).is_err());
        assert!(OutputTimes::lin_log(2.0, 1.0, 4, 3).is_err());

        assert_eq!(extension_schedule(None, 0, 2.0), Ok(vec![2.0]));
        assert_eq!(extension_schedule(Some(0.5), 2, 3.0), Ok(vec![0.5, 0.5, 2.0]));
        assert!(extension_schedule(None, 2, 3.0).is_err());
        assert!(extension_schedule(Some(0.5), 0, 3.0).is_err());
        assert!(extension_schedule(Some(2.0), 2, 3.0).is_err());
    }

    #[test]
    fn sample_at_every_output_time() {
        let emodel = Arc::new(Fuzzyfold::default());
        let rmodel = Arrhenius::new(emodel.temperature(), 1e5);
        let sequence = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let start = MultiPairTable::try_from("..........").unwrap();
        let num_ext = sequence.len() - start.len();
        let moves = LoopNeighbors::try_from((sequence.clone(), &start, emodel, MoveTypes::NONE)).unwrap();

        let out = OutputTimes::transcription(1e-3, num_ext, 1e-2).unwrap();
        let schedule = extension_schedule(Some(1e-3), num_ext, out.end()).unwrap();
        let mut ssa = SSA::from((moves, rmodel));
        let lengths = ssa.sample(&mut StdRng::seed_from_u64(7), &schedule, &out,
            |m| m.current_structure_length());
        assert_eq!(lengths.len(), out.times().len());
        assert_eq!(lengths[0], start.len());
        assert_eq!(*lengths.last().unwrap(), sequence.len());
        assert!(lengths.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn test_flux_after_moves() {
        let emodel = Arc::new(Fuzzyfold::default());
        let rmodel = Arrhenius::new(emodel.temperature(), 1.0).with_three_way(1.0).with_four_way(1.0);
        let policy = MoveTypes::SHIFTS;

        let sequence = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").expect("Invalid sequence?");
        let pairings =     MultiPairTable::try_from(".((((....)))).((((........))))...............").expect("Invalid structure?");

        let ltab = LoopTable::try_from((sequence.clone(), &pairings, emodel.clone()))
            .expect("Invalid sequence/structure combination");
        let moves = LoopNeighbors::from((ltab, policy));
 
        let mut simulator = SSA::from((moves, rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng(), 100., |_, _, _, w| {
            steps += 1;
            let moves1 = w.propose_moves().collect::<Vec<_>>();
            let p = MultiPairTable::try_from(&w.current_structure()).expect("Invalid structure?");
            let ltab = LoopTable::try_from((sequence.clone(), &p, emodel.clone()))
                .expect("Invalid sequence/structure combination");
            let moves = LoopNeighbors::from((ltab, policy));
            let moves2 = moves.propose_moves().collect::<Vec<_>>();
            assert!(same_moves(&moves1, &moves2));
            true
        });
        assert!(steps > 0, "Simulation must perform at least one step");
    }


    #[test]
    fn test_flux_after_moves_circular() {
        // On circular strands, the moves after each step are the moves of a
        // table built from scratch, and the energy is that of the structure.
        let emodel = Arc::new(Fuzzyfold::default());
        let rmodel = Arrhenius::new(emodel.temperature(), 1.0).with_three_way(1.0).with_four_way(1.0);
        let policy = MoveTypes::SHIFTS;
        for (seq, db) in [
            ("GGGAAACCCA^", "..........^"),
            ("GCGCAAAAGCGCAAAAAGCGCUUUUGCGC^", ".............................^"),
            ("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC^", 
             ".((((....)))).((((........))))...............^"),
        ] {
            let sequence = NucleotideVec::try_from(seq).unwrap();
            let pairings = MultiPairTable::try_from(db).unwrap();
            let ltab = LoopTable::try_from((sequence.clone(), &pairings, emodel.clone())).unwrap();
            let mut simulator = SSA::from((LoopNeighbors::from((ltab, policy)), rmodel));
            let mut steps = 0;
            simulator.simulate(&mut StdRng::seed_from_u64(7), 1e9, |_, _, _, w| {
                steps += 1;
                let structure = w.current_structure();
                let p = MultiPairTable::try_from(&structure).unwrap();
                assert_eq!(w.current_energy(), emodel.energy_of_structure(&sequence, &p).unwrap());
                let ltab = LoopTable::try_from((sequence.clone(), &p, emodel.clone())).unwrap();
                let moves = LoopNeighbors::from((ltab, policy)).propose_moves().collect::<Vec<_>>();
                assert!(same_moves(&w.propose_moves().collect::<Vec<_>>(), &moves), "{structure}");
                steps < 2000
            });
            assert!(steps > 100, "{seq}: {steps} steps");
        }
    }

    #[test]
    fn test_flux_after_moves_cotranscriptional() {
        // After each step and each extension, the moves are those of a table
        // built from scratch, with and without shift moves.
        let emodel = Arc::new(Fuzzyfold::default());
        let seq = "UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC";
        let sequence = NucleotideVec::try_from(seq).unwrap();
        let start = MultiPairTable::try_from("..........").unwrap();
        for shifts in [MoveTypes::NONE, MoveTypes::SHIFTS] {
            let k = if shifts == MoveTypes::SHIFTS { 1.0 } else { 0.0 };
            let rmodel = Arrhenius::new(emodel.temperature(), 1.0).with_three_way(k).with_four_way(k);
            let moves = LoopNeighbors::try_from((sequence.clone(), &start, emodel.clone(), shifts)).unwrap();
            let mut simulator = SSA::from((moves, rmodel));
            let times = vec![5.0; seq.len() - start.len() + 1];
            let mut steps = 0;
            simulator.co_simulate(&mut StdRng::seed_from_u64(3), &times, |_, _, _, w| {
                steps += 1;
                let structure = w.current_structure();
                let p = MultiPairTable::try_from(&structure).unwrap();
                let fresh = LoopNeighbors::try_from((sequence.clone(), &p, emodel.clone(), shifts)).unwrap();
                let moves = fresh.propose_moves().collect::<Vec<_>>();
                assert!(same_moves(&w.propose_moves().collect::<Vec<_>>(), &moves), "{structure}");
                true
            });
            assert!(steps > 100, "{shifts:?}: {steps} steps");
        }
    }

    #[test]
    fn test_ends_in_the_last_reported_structure() {
        // The waiting time belongs to the structure the callback sees: when
        // it ends after t_max, the structure stays; also at each extension.
        let emodel = Arc::new(Fuzzyfold::default());
        let rmodel = Arrhenius::new(emodel.temperature(), 1.0).with_three_way(1.0).with_four_way(1.0);
        let seq = "UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC";
        let sequence = NucleotideVec::try_from(seq).unwrap();
        for (seed, t_max) in [(1, 0.5), (2, 3.0), (3, 20.0)] {
            let start = MultiPairTable::try_from(".".repeat(seq.len()).as_str()).unwrap();
            let moves = LoopNeighbors::try_from((sequence.clone(), &start, emodel.clone(), MoveTypes::SHIFTS)).unwrap();
            let mut simulator = SSA::from((moves, rmodel));
            let mut last = None;
            simulator.simulate(&mut StdRng::seed_from_u64(seed), t_max, |_, _, _, w| {
                last = Some(w.current_structure());
                true
            });
            assert_eq!(Some(simulator.current_structure()), last, "t_max {t_max}");
        }
        // Before each extension: the structure of the last callback, plus an
        // unpaired base.
        let start = MultiPairTable::try_from("..........").unwrap();
        let moves = LoopNeighbors::try_from((sequence.clone(), &start, emodel.clone(), MoveTypes::SHIFTS)).unwrap();
        let mut simulator = SSA::from((moves, rmodel));
        let mut previous: Option<DotBracketVec> = None;
        let times = vec![2.0; seq.len() - start.len() + 1];
        simulator.co_simulate(&mut StdRng::seed_from_u64(4), &times, |_, _, _, w| {
            let structure = w.current_structure();
            if let Some(p) = &previous && p.len() < structure.len() {
                assert_eq!(&structure.to_string()[..p.len()], p.to_string(), "{structure}");
                assert!(structure.to_string().ends_with('.'));
            }
            previous = Some(structure);
            true
        });
    }

    #[test]
    fn test_simple_ssa_simulation() {
        setup_ssa_input!(moves, rmodel, "GUACACGUCG", "..........");
        let mut rng = StdRng::seed_from_u64(42);
        let mut simulator = SSA::from((moves, rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng, 1.0, |_, _, _, _| {
            steps += 1;
            true
        });
        assert!(steps > 0, "Simulation must perform at least one step");
    }

    #[test]
    fn test_cotr_ssa_simulation() {
        setup_ssa_input!(moves, rmodel, "GUACACGUCG", "......");
        let mut rng = StdRng::seed_from_u64(42);
        let mut simulator = SSA::from((moves, rmodel));

        let mut last = 0.0;
        simulator.co_simulate(&mut rng, &[4000.0; 5], |t, tinc, _, _| {
            assert!(t >= last && tinc >= 0.0);
            last = t;
            true
        });
        assert_eq!(simulator.move_set().current_structure_length(), 10);
    }
}



