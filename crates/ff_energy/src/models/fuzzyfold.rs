//! The energy model of fuzzyfold: the loops of ViennaRNA, and pseudoknots
//! with the pseudoloop terms of `PKParams` and the scaled loops of strained
//! stems.

use std::fmt;
use std::sync::atomic::{AtomicBool, Ordering};
use log::warn;

use crate::Base;
use crate::EnergyModel;
use crate::EnergyError;
use crate::parameters::*;
use crate::Pseudoknot;
use crate::PseudoknotLoop;
use crate::Gap;
use crate::NestedLoop;
use crate::Ring;
use crate::RingDecomposition;
use crate::circular_strands;
use super::viennarna::ViennaRNA;

/// The energy model of fuzzyfold: the loops of [`ViennaRNA`], pseudoknots
/// with `PKParams`, and structures of several strands independent of the
/// order of the strands. The open circle costs the initiation of a hairpin
/// loop of its ring (0 in ViennaRNA). This is the model under development; 
/// `ViennaRNA` stays as it is.
pub struct Fuzzyfold {
    /// The loop energies, shared with ViennaRNA until the models diverge.
    vrna: ViennaRNA,
    /// The pseudoknot parameters.
    pk_params: PKParams,
}

impl Default for Fuzzyfold {
    fn default() -> Self {
        Self::from_thermo_params(&RNA_EXTENDED, 37.0)
    }
}

/// The strands joined into complexes, counted while the rings are
/// evaluated. The complexes are the components of the loop graph (a vertex
/// per loop, an edge per pair): loops minus pairs plus the cycles of the
/// pseudoknots. Each linear strand starts one open loop; circular strands
/// are counted from the sequence.
struct StrandJoins {
    strands: usize,
    loops: usize,
    pairs: usize,
    cycles: usize,
}

impl StrandJoins {
    fn new(circular: usize) -> Self {
        Self { strands: circular, loops: 0, pairs: 0, cycles: 0 }
    }

    fn count(&mut self, l: &NestedLoop) {
        self.strands += matches!(l.gap(), Gap::Open { .. }) as usize;
        self.loops += 1;
        self.pairs += l.pairs().len();
    }

    fn add(&mut self, ring: &Ring) {
        match ring {
            Ring::Loop(l) => self.count(l),
            Ring::Pseudoknot(pk) => {
                pk.loops.iter().for_each(|l| self.count(&l.nested));
                pk.stacks.iter().chain(&pk.interior_loops).for_each(|l| self.count(l));
                self.cycles += pk.cycles;
            }
        }
    }

    /// Strands minus complexes; each pair is seen from two loops.
    fn joins(&self) -> i32 {
        if self.strands < 2 {
            return 0;
        }
        (self.strands + self.pairs / 2 - self.loops - self.cycles) as i32
    }
}

impl Fuzzyfold {
    /// A model from fitted parameters, at 37 C only.
    pub fn from_andronescu_params(params: &'static AndronescuParams) -> Self {
        Self { vrna: ViennaRNA::from_andronescu_params(params), pk_params: PKParams::default() }
    }

    /// A model from thermodynamic parameters, rescaled to a temperature in
    /// Celsius. The pseudoknot parameters are not rescaled.
    pub fn from_thermo_params(params: &'static RNAThermoParams, celsius: f64) -> Self {
        Self { vrna: ViennaRNA::from_thermo_params(params, celsius), pk_params: PKParams::default() }
    }

    /// Replace the pseudoknot parameters (see [`PKParams`] for how they are
    /// applied). Every model starts with `PKParams::default()`.
    ///
    /// The parameters are given at 37 C only: they are not rescaled with
    /// the temperature of the model.
    ///
    /// ```rust
    /// use ff_energy::{Fuzzyfold, parameters::{RNA_TURNER_2004, PKParams}};
    /// let model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0)
    ///     .with_pseudoknot_params(PKParams { stem: 400, ..PKParams::default() });
    /// ```
    pub fn with_pseudoknot_params(mut self, params: PKParams) -> Self {
        self.pk_params = params;
        self
    }

    /// The energy of each loop and pseudoknot of a structure, and of the
    /// duplex initiations of its complexes (level 0), and for a pseudoknot
    /// also its terms (levels 1 and 2):
    /// `(level, description, energy)`, in the order of the decomposition.
    /// The energies at level 0 sum to the energy of the structure, the terms
    /// of a pseudoknot to its energy.
    pub fn energy_terms<T: RingDecomposition>(&self, sequence: &[Base], structure: &T
    ) -> Result<Vec<(usize, String, i32)>, EnergyError> {
        // The terms of a pseudoknot come before its total: list them after.
        let (mut terms, mut inner) = (Vec::new(), Vec::new());
        self.sum_of_rings(sequence, structure, &mut |level, args, e| {
            let term = (level, args.to_string(), e);
            if level == 0 {
                terms.push(term);
                terms.append(&mut inner);
            } else {
                inner.push(term);
            }
        })?;
        Ok(terms)
    }

    /// The sum of the rings, and the duplex initiations of the complexes.
    /// Each ring and term is also passed to `note`, with its level.
    fn sum_of_rings<T: RingDecomposition, N: FnMut(usize, fmt::Arguments, i32)>(
        &self, sequence: &[Base], structure: &T, note: &mut N,
    ) -> Result<i32, EnergyError> {
        let (mut total, mut joins) = (Ok(0), StrandJoins::new(circular_strands(sequence)));
        structure.for_each_ring(|r| {
            let Ok(sum) = total else { return };
            joins.add(r);
            let en = match r {
                Ring::Loop(l) => self.energy_of_nested_loop(sequence, l),
                Ring::Pseudoknot(pk) => self.pseudoknot(sequence, pk, note),
            };
            total = en.map(|en| {
                note(0, format_args!("{r}"), en);
                sum + en
            });
        });
        let joins = joins.joins();
        if total.is_ok() && joins > 0 {
            note(0, format_args!("Duplex initiation x {joins}"), self.vrna.duplex_init * joins);
        }
        total.map(|en| en + self.vrna.duplex_init * joins)
    }

    /// The energy of a `Pseudoknot`: initiation per knot, Pb per strained
    /// stem, the scaled interior loops of the strained stems, Pup per unpaired
    /// base of a strained segment, Pps per branch in a strained segment, and
    /// the helix ends of its loops: terminal penalty, and a mismatch if both
    /// flanks are in the same context (`mismatch_multi` if strained,
    /// `mismatch_exterior` if not), else the dangle on each side. The duplex
    /// initiation is a term of the complex (see `energy_of_structure`). Each
    /// term is also passed to `note`, with its level.
    fn pseudoknot(&self, sequence: &[Base], pk: &Pseudoknot,
        note: &mut dyn FnMut(usize, fmt::Arguments, i32),
    ) -> Result<i32, EnergyError> {
        let (p, vrna) = (&self.pk_params, &self.vrna);
        if (vrna.temperature - T_REF_CELSIUS).abs() > 1e-6 {
            warn_unscaled_pseudoknots(vrna.temperature);
        }
        let mut en = 0;
        let mut add = |level: usize, args: fmt::Arguments, e: i32| {
            note(level, args, e);
            en += e;
        };
        add(1, format_args!("initiation x {}", pk.knots), p.init * pk.knots as i32);
        add(1, format_args!("strained stems x {}", pk.stems), p.stem * pk.stems as i32);
        let scaled = pk.stacks.iter().map(|l| (l, p.stack_scale))
            .chain(pk.interior_loops.iter().map(|l| (l, p.interior_scale)));
        for (l, scale) in scaled {
            let e = (scale * vrna.nested_loop(sequence, l)? as f64).round() as i32;
            // Named here, not by the Display of the loop, which decides by
            // index: the backbone may pass the wraps of circles.
            let [(a, b), (c, d)] = l.pairs() else { unreachable!("A loop of a stem has two pairs.") };
            add(2, format_args!("{:<8} ({:>3}, {:>3}), ({:>3}, {:>3}) x{scale}", "Interior",
                c.min(d), c.max(d), a.min(b), a.max(b)), e);
        }
        for PseudoknotLoop { nested, strained, stem_ends } in &pk.loops {
            let pairs = nested.pairs();
            let segments: Vec<_> = nested.segments(sequence).collect();
            let open = matches!(nested.gap(), Gap::Open { .. });
            debug_assert!(!open || !strained[0] && !strained[pairs.len()],
                "The segments to the strand ends are not strained.");
            // The bases of the strained segments, without their paired ends.
            let unpaired: usize = segments.iter().zip(strained)
                .filter(|&(_, &s)| s)
                .map(|(s, _)| s.len() - 2)
                .sum();
            add(1, format_args!("{nested}, unpaired x {unpaired}"), p.unpaired * unpaired as i32);
            for (t, &(a, b)) in pairs.iter().enumerate() {
                let u = (t + 1) % segments.len();
                // The flanks of the pair in the loop: the base before it in
                // the segment before, the base after it in the segment after.
                let (before, after) = (&segments[t], &segments[u]);
                let d5 = before.len().checked_sub(2).map(|d| before[d]);
                let d3 = after.get(1).copied();
                let pair = (sequence[a], sequence[b]);
                let (context, e) = match (strained[t], strained[u]) {
                    (true, true) if stem_ends[t] =>
                        ("strained", vrna.helix_end(pair, d5, d3, &vrna.mismatch_multi)?),
                    (true, true) =>
                        ("branch", p.branch + vrna.helix_end(pair, d5, d3, &vrna.mismatch_multi)?),
                    (false, false) =>
                        ("exterior", vrna.helix_end(pair, d5, d3, &vrna.mismatch_exterior)?),
                    _ => {
                        debug_assert!(stem_ends[t], "Only a strained stem ends a strained segment.");
                        ("mixed", vrna.helix_end(pair, d5, None, &vrna.mismatch_exterior)?
                            + vrna.dangles(pair, None, d3, &vrna.mismatch_exterior))
                    }
                };
                add(2, format_args!("helix end ({a}, {b}), {context}"), e);
            }
        }
        Ok(en)
    }

    /// The open circle closes a ring of n bases: the initiation of a hairpin
    /// loop of the same ring, n - 2 unpaired bases closed by a pair.
    fn open_circle(&self, n: usize) -> Result<i32, EnergyError> {
        self.vrna.hairpin_initiation(n.saturating_sub(2))
    }
}

impl EnergyModel for Fuzzyfold {
    fn temperature(&self) -> f64 {
        self.vrna.temperature()
    }

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        self.vrna.can_pair(b1, b2)
    }

    fn min_hairpin_size(&self) -> usize {
        self.vrna.min_hairpin_size()
    }

    /// The sum of the rings, and one duplex initiation for each strand that
    /// joins a complex (strands minus complexes): a strand may be joined by
    /// a pseudoknot alone, without a loop that joins two strands.
    fn energy_of_structure<T: RingDecomposition>(&self, sequence: &[Base], structure: &T
    ) -> Result<i32, EnergyError> {
        self.sum_of_rings(sequence, structure, &mut |_, _, _| ())
    }

    /// Without the duplex initiation of a loop that joins two strands: it is
    /// a term of the complex (see `energy_of_structure`).
    fn energy_of_nested_loop(&self, sequence: &[Base], nested: &NestedLoop
    ) -> Result<i32, EnergyError> {
        match nested.gap() {
            Gap::Wrap { last, first } if nested.pairs().is_empty() => self.open_circle(last - first + 1),
            _ => self.vrna.nested_loop(sequence, nested),
        }
    }

    fn energy_of_pseudoknot(&self, sequence: &[Base], pseudoknot: &Pseudoknot
    ) -> Result<i32, EnergyError> {
        self.pseudoknot(sequence, pseudoknot, &mut |_, _, _| ())
    }

}

/// Pseudoknot parameters are given at 37 C only. Warn (once per process) 
/// when a pseudoknot is evaluated at another temperature.
fn warn_unscaled_pseudoknots(celsius: f64) {
    static WARNED: AtomicBool = AtomicBool::new(false);
    if !WARNED.swap(true, Ordering::Relaxed) {
        warn!("Pseudoknot parameters are given at 37 C and are not rescaled to {} C.", celsius);
    }
}

#[cfg(test)]
mod tests {
    use ff_structure::MultiPairTable;
    use crate::NucleotideVec;
    use super::*;

    /// The energy of a structure.
    fn energy(model: &Fuzzyfold, seq: &str, db: &str) -> Result<i32, EnergyError> {
        model.energy_of_structure(&NucleotideVec::try_from_rna(seq).unwrap(), &MultiPairTable::try_from(db).unwrap())
    }

    /// The pseudoknot of a structure.
    fn pseudoknot_of(db: &str) -> Pseudoknot {
        MultiPairTable::try_from(db).unwrap().rings().into_iter()
            .find_map(|r| match r { Ring::Pseudoknot(pk) => Some(pk), _ => None }).unwrap()
    }

    /// Terminal penalty and mismatch of a helix end (p5, p3) with both 
    /// flanks in the pseudoloop.
    fn multi(model: &Fuzzyfold, seq: &[Base], (p5, p3): (usize, usize)) -> i32 {
        model.vrna.helix_end((seq[p5], seq[p3]), Some(seq[p5 - 1]), Some(seq[p3 + 1]), 
            &model.vrna.mismatch_multi).unwrap()
    }

    /// Terminal penalty and the dangles d5 and d3 (if any) of a helix end 
    /// (p5, p3) with flanks in different contexts.
    fn mixed(model: &Fuzzyfold, seq: &[Base], (p5, p3): (usize, usize), d5: Option<usize>, d3: Option<usize>) -> i32 {
        let pair = (seq[p5], seq[p3]);
        let m = &model.vrna.mismatch_exterior;
        model.vrna.helix_end(pair, d5.map(|d| seq[d]), None, m).unwrap() 
            + model.vrna.dangles(pair, None, d3.map(|d| seq[d]), m)
    }

    #[test]
    fn pseudoknot_with_closing_pair() {
        // (.[[[....(((((]]]......)))))..): the knot's outside is a stretch of 
        // its pseudoloop, the closing pair a child; initiated as in the exterior.
        let seq = NucleotideVec::try_from_rna("GAGCGAUUUCUGACCGCUUUUUUGUCAGAAC").unwrap();
        let pk_loop = pseudoknot_of("(.[[[....(((((]]]......)))))..)");
        assert_eq!((pk_loop.knots, pk_loop.stems), (1, 2));
        assert!(pk_loop.loops.iter().all(|l| l.strained.iter().all(|&s| s)));

        let model = Fuzzyfold::default();
        let pk = PKParams::default();
        // Each stack is scaled and rounded on its own.
        let stack = |(i, j): (usize, usize), (p, q): (usize, usize)| (pk.stack_scale
            * model.vrna.interior(&seq[i..=p], &seq[q..=j]).unwrap() as f64).round() as i32;
        let stacks = stack((2, 16), (3, 15)) + stack((3, 15), (4, 14))
            + (9..13).map(|k| stack((k, 36 - k), (k + 1, 35 - k))).sum::<i32>();
        // Helix ends, all in the pseudoloop: two per stem, and the closing pair.
        let ends = [(2, 16), (14, 4), (9, 27), (23, 13), (30, 0)].iter()
            .map(|&e| multi(&model, &seq, e)).sum::<i32>();
        // Unpaired: 1 | 5..=8 | - | 17..=22 | 28, 29
        let expected = pk.init + pk.branch + 2 * pk.stem + 13 * pk.unpaired + stacks + ends;
        assert_eq!(model.energy_of_pseudoknot(&seq, &pk_loop).unwrap(), expected);
    }

    #[test]
    fn pseudoknot_h_type() {
        // ((([[[)))...]]]
        let seq = NucleotideVec::try_from_rna("GGGCCCCCCAAAGGG").unwrap();
        let pk_loop = pseudoknot_of("((([[[)))...]]]");
        let model = Fuzzyfold::default();
        let pk = PKParams::default();
        // Each stack is scaled and rounded on its own.
        let stack = |seq: &[Base], (i, j): (usize, usize), (p, q): (usize, usize)| (pk.stack_scale
            * model.vrna.interior(&seq[i..=p], &seq[q..=j]).unwrap() as f64).round() as i32;
        let stacks = stack(&seq, (0, 8), (1, 7)) + stack(&seq, (1, 7), (2, 6)) 
                   + stack(&seq, (3, 14), (4, 13)) + stack(&seq, (4, 13), (5, 12));
        // The outer ends have one flank at an end of the sequence.
        let ends = mixed(&model, &seq, (0, 8), None, Some(9)) + multi(&model, &seq, (6, 2))
            + mixed(&model, &seq, (3, 14), Some(2), None) + multi(&model, &seq, (12, 5));
        let expected = pk.init + 2 * pk.stem + 3 * pk.unpaired + stacks + ends;
        assert_eq!(model.energy_of_pseudoknot(&seq, &pk_loop).unwrap(), expected);

        // A pseudoknot in a stretch of another one is part of its pseudoloop:
        // ((([[[)))..<<<{{{>>>}}}..]]] has four stems, no children, and one 
        // initiation. (Counting crossing components instead depends on the 
        // drawing, see pseudoknot_outside_is_a_stretch.)
        let seq = NucleotideVec::try_from_rna("GGGCCCCCCAAGGGCCCCCCGGGAAGGG").unwrap();
        let pk_loop = &pseudoknot_of("((([[[)))..<<<{{{>>>}}}..]]]");
        assert_eq!((pk_loop.knots, pk_loop.stems), (1, 4));
        let stacks: i32 = [(0, 8), (3, 27), (11, 19), (14, 22)].iter()
            .flat_map(|&(i, j)| (0..2).map(move |t| ((i + t, j - t), (i + t + 1, j - t - 1))))
            .map(|(c, i)| stack(&seq, c, i))
            .sum();
        let ends = mixed(&model, &seq, (0, 8), None, Some(9)) + multi(&model, &seq, (6, 2))
            + mixed(&model, &seq, (3, 27), Some(2), None) + multi(&model, &seq, (25, 5))
            + [(11, 19), (17, 13), (14, 22), (20, 16)].iter().map(|&e| multi(&model, &seq, e)).sum::<i32>();
        // Unpaired: 9, 10 | 23, 24
        let expected = pk.init + 4 * pk.stem + 4 * pk.unpaired + stacks + ends;
        assert_eq!(model.energy_of_pseudoknot(&seq, pk_loop).unwrap(), expected);
    }

    #[test]
    fn pseudoknot_cut_invariance() {
        // Two strands, two parallel helices X and Y: a genuine H-type in both 
        // cuts. The overhangs end up in the outside in one cut and in the 
        // nicked stretch in the other; the energy must not change.
        //   A = AAAA GGGC AAA CGGG AA,  B = A GCCC UUU CCCG
        let model = Fuzzyfold::default();
        let energy = |seq: &str, db: &str| energy(&model, seq, db).unwrap();
        assert_eq!(
            energy("AAAAGGGCAAACGGGAA+AGCCCUUUCCCG", "....((((...[[[[..+.))))...]]]]"),
            energy("AGCCCUUUCCCG+AAAAGGGCAAACGGGAA", ".((((...[[[[+....))))...]]]].."),
        );
    }

    #[test]
    fn pseudoknot_strand_order_invariance() {
        // Three strands, three helices: X = A5'.B5', Y = A3'.C5', Z = B3'.C3'.
        //   A = GGGC AAA CGGG,  B = GCCC AAA GCGC,  C = CCCG AAA GCGC
        // Y ends blunt at the A|C junction, so it is free in every drawing; 
        // what remains is the H-type X x Z.
        let model = Fuzzyfold::default();
        let energy = |seq: &str, db: &str| energy(&model, seq, db).unwrap();
        // Kissing hairpin, a band split by a nick, and the H-type with Y nested.
        let abc = energy("GGGCAAACGGG+GCCCAAAGCGC+CCCGAAAGCGC", "((((...[[[[+))))...((((+]]]]...))))");
        let cab = energy("CCCGAAAGCGC+GGGCAAACGGG+GCCCAAAGCGC", "[[[[...((((+((((...]]]]+))))...))))");
        let acb = energy("GGGCAAACGGG+CCCGAAAGCGC+GCCCAAAGCGC", "((((...((((+))))...[[[[+))))...]]]]");
        assert_eq!(abc, acb);
        assert_eq!(cab, acb);

        // Nested in one order, pseudoknotted in the others: a multiloop closed 
        // by B1-C6 with two branches, each of which hides a nick.
        let nested = energy("GUUGGG+CCCAAUAGGACAA+CCGGGUAU", ".(((((+))))...((....+))....).");
        assert_eq!(energy("CCCAAUAGGACAA+GUUGGG+CCGGGUAU", "((((...[[....+.{))))+]]....}."), nested);
        assert_eq!(energy("CCGGGUAU+CCCAAUAGGACAA+GUUGGG", "((....[.+{{{{...))....+.]}}}}"), nested);

        // A 3-way junction between H1, H2 and H3 is strained in every
        // drawing: H1 and H2 are two stems either way.
        let abc = energy("CUGUCUCUG+UUGAUUACACAGU+AACAGACGGCGCACAC", "((((.[[[[+.{{.....))))(+).]]]].......}}.");
        let acb = energy("CUGUCUCUG+AACAGACGGCGCACAC+UUGAUUACACAGU", "((((.((((+[.)))).......((.+.)).....))))]");
        assert_eq!(abc, acb);
    }

    #[test]
    fn pseudoknot_outside_is_a_stretch() {
        // The outside of a knot is just another stretch of its pseudoloop, 
        // whichever face the drawing makes the outside.
        let model = Fuzzyfold::default();
        let energy = |seq: &str, db: &str| energy(&model, seq, db).unwrap();

        // The knot in the exterior, or enclosed by a multiloop (absorbed).
        assert_eq!(
            energy("GACAAGGCUAUGGAA+AUUGAAAUUCCUUG", ".([[[....((....+......)))..]]]"),
            energy("AUUGAAAUUCCUUG+GACAAGGCUAUGGAA", "......(((..[[[+.)]]]....))...."),
        );
        // One knot with four stems, or a knot in a stretch of another (merged).
        assert_eq!(
            energy("AGACGACUCACCUUU+ACACAAGG+UGGACUGGAACACCAA+GUUACUGAA", 
                   "((.........(((.+[[..))).+....)).......((.+]]...)).."),
            energy("ACACAAGG+AGACGACUCACCUUU+UGGACUGGAACACCAA+GUUACUGAA", 
                   "((..(((.+[[.........))).+....]].......[[.+))...]].."),
        );
    }

    #[test]
    fn circular_rotation_invariance() {
        // The strand wrap may lie in any loop: every rotation of a circle has
        // the same energy. Three hairpins: the exterior loop is a multiloop.
        let seq = "AAUGCUAAGACAUUUCCCUUCAGGGGGGGCUCCCCCGC";
        let db  = "((((......)))).(((....)))((((...))))..";
        let n = seq.len();
        let mut partner = vec![None; n];
        let mut stack = Vec::new();
        for (k, c) in db.chars().enumerate() {
            match c {
                '(' => stack.push(k),
                ')' => { let i = stack.pop().unwrap(); partner[i] = Some(k); partner[k] = Some(i); }
                _ => (),
            }
        }
        let model = Fuzzyfold::default();
        let energy = |r: usize| {
            let rseq: String = seq[r..].chars().chain(seq[..r].chars()).chain(['^']).collect();
            let rdb: String = (0..n).map(|k| match partner[(k + r) % n] {
                Some(p) if (p + n - r) % n > k => '(',
                Some(_) => ')',
                None => '.',
            }).chain(['^']).collect();
            model.energy_of_structure(&NucleotideVec::try_from_rna(&rseq).unwrap(),
                &MultiPairTable::try_from(rdb.as_str()).unwrap()).unwrap()
        };
        let e0 = energy(0);
        for r in 1..n {
            assert_eq!(energy(r), e0, "rotation {r}");
        }
    }

    #[test]
    fn open_circle_energy() {
        // The initiation of a hairpin loop with n - 2 unpaired bases (Turner
        // 2004: 7.10 kcal/mol for 18, 7.70 + 1.07856 ln(98 / 30) for 98).
        let model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let energy = |n: usize| {
            let seq = format!("{}^", "A".repeat(n));
            let db = format!("{}^", ".".repeat(n));
            model.energy_of_structure(&NucleotideVec::try_from_rna(&seq).unwrap(),
                &MultiPairTable::try_from(db.as_str()).unwrap())
        };
        assert_eq!(energy(20).unwrap(), 710);
        assert_eq!(energy(100).unwrap(), 897);
        assert_eq!(energy(5).unwrap(), 540);
        assert!(matches!(energy(4), Err(EnergyError::HairpinTooSmall { size: 2, min: 3 })));
        // ViennaRNA keeps the open circle at 0.
        let vrna = ViennaRNA::default();
        assert_eq!(vrna.energy_of_structure(&NucleotideVec::try_from_rna("AAAAAAAAAA^").unwrap(),
            &MultiPairTable::try_from("..........^").unwrap()).unwrap(), 0);
    }

    #[test]
    fn circular_pseudoknots() {
        let model = Fuzzyfold::default();
        let energy = |seq: &str, db: &str| energy(&model, seq, db);
        // An H-type in the loop through the wrap, cut anywhere.
        let seq = "GCGAUUUCUGACCGCUUUUUUGUCAG";
        let db = "[[[....(((((]]]......)))))";
        let e0 = energy(&format!("{seq}^"), &format!("{db}^")).unwrap();
        assert_eq!(energy(&format!("{}{}^", &seq[1..], &seq[..1]), "((....[[[[[))(......]]]]])^").unwrap(), e0);
        // A pseudoknot enclosed by a pair: the circle is a hairpin around it.
        let seq = "AAGGAGCGAUUUCUGACCGCUUUUUUGUCAGAACCAA";
        let db  = "..((.[[[....(((((]]]......)))))..))..";
        let linear = energy(seq, db).unwrap();
        let circular = energy(&format!("{seq}^"), &format!("{db}^")).unwrap();
        let bases = NucleotideVec::try_from_rna(seq).unwrap();
        let exterior = model.vrna.exterior([&bases[0..=2], &bases[34..=36]]).unwrap();
        // The loop through the wrap: C34 A35 A36 | A0 A1 G2.
        let hairpin = model.vrna.hairpin(&NucleotideVec::try_from_rna("CAAAAG").unwrap()).unwrap();
        assert_eq!(circular, linear - exterior + hairpin);
    }
}
