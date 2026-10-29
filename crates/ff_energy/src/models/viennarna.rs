use paste::paste;

use crate::Base;
use crate::PairTypeRNA;
use crate::EnergyModel;
use crate::EnergyError;
use crate::parameters::*;
use crate::Ring;
use crate::RingDecomposition;
use crate::Pseudoknot;
use crate::Gap;
use crate::NestedLoop;

/// The energy model of ViennaRNA 2.6 with its standard settings (dangles on
/// both sides, special hairpins, `RNAeval -d2`), for any of the parameter
/// sets in [`crate::parameters`]. Pseudoknots are not supported. Its loop
/// energies are shared with [`crate::Fuzzyfold`] until the models diverge.
///
/// At temperatures other than 37 C, the rescaled parameters are rounded;
/// ViennaRNA truncates them, so energies differ by a few tenths of a
/// kcal/mol there.
pub struct ViennaRNA {
    /// Steric constraint on minimum hairpin length. (Always three.)
    min_hp_size: usize,

    /// Intramolecular initiation free energy.
    pub(super) duplex_init: i32,

    /// The temperature of measured parameters. 
    /// (Nonensical for fitted parameters.)
    pub(super) temperature: f64, 

    /// An extended parameter table for stacks. 
    stack: ExtendedStackParams,

    mismatch_hairpin: MismatchParams,
    mismatch_interior: MismatchParams,
    mismatch_interior_1n: MismatchParams,
    mismatch_interior_23: MismatchParams,
    pub(super) mismatch_multi: MismatchParams,
    pub(super) mismatch_exterior: MismatchParams,
    dangle5: DangleParams,
    dangle3: DangleParams,
    int11: Int11Params,
    int21: Int21Params,
    int22: Int22Params,

    /// Paramters for loops with len <= 30.
    hairpin: LoopParams,
    bulge: LoopParams,
    interior: LoopParams,

    /// Extrapolation constant for loops with len > 30 based on polymer theory.
    lxc: f64,

    /// An extended parameter table for closing penalties. 
    closing_penalty: ClosingParams,

    /// Asymmetric internal loop correction.
    ninio: i32,
    ninio_max: i32,

    /// Multiloop penalty for unpaired bases.
    ml_base: i32,
    /// Multiloop initiation penalty.
    ml_closing: i32,
    /// Multiloop penalty per outgoing stem.
    ml_intern: i32,

    triloops: Vec<LoopEntry>,
    tetraloops: Vec<LoopEntry>,
    hexaloops: Vec<LoopEntry>,

}

macro_rules! rescale_params {
    ($field:ident, $params:ident, $scale:ident) => {
        paste! {
                $params.[<$field _en37>]
                    .rescale($params.[<$field _enth>], $scale)
        }
    };
}

macro_rules! rescale_param {
    ($field:ident, $params:ident, $scale:ident) => {
        paste! {
                $params.[<$field _en37>]
                    .rescale(&$params.[<$field _enth>], $scale)
        }
    };
}


impl Default for ViennaRNA {
    fn default() -> Self {
        Self::from_thermo_params(&RNA_TURNER_2004, 37.0)
    }
}

impl ViennaRNA {
    /// Initializes a model from fitted parameters, which means there
    /// is no possiblity to change the temperature.
    pub fn from_andronescu_params(params: &'static AndronescuParams) -> Self {
        Self {
            min_hp_size: 3,
            temperature: 37.0,

            stack: *params.stack,
            mismatch_hairpin: *params.mismatch_hairpin,
            mismatch_interior: *params.mismatch_interior,
            mismatch_interior_1n: *params.mismatch_interior_1n,
            mismatch_interior_23: *params.mismatch_interior_23,
            mismatch_multi: *params.mismatch_multi,
            mismatch_exterior: *params.mismatch_exterior,

            dangle5: *params.dangle5,
            dangle3: *params.dangle3,

            int11: *params.int11,
            int21: *params.int21,
            int22: *params.int22,

            hairpin: *params.hairpin,
            bulge: *params.bulge,
            interior: *params.interior,

            duplex_init: params.duplex_init,
            closing_penalty: *params.closing,

            lxc: params.lxc,

            ninio: params.ninio,
            ninio_max: params.ninio_max,

            ml_base: params.ml_base,
            ml_closing: params.ml_closing,
            ml_intern: params.ml_intern,

            triloops: params.triloops.to_vec(),
            tetraloops: params.tetraloops.to_vec(),
            hexaloops: params.hexaloops.to_vec(),
        }
    }

    /// Initializes a model from thermodynamic parameters. That's the default!
    pub fn from_thermo_params(params: &'static RNAThermoParams, celsius: f64) -> Self {
        if (celsius - T_REF_CELSIUS).abs() < 1e-6 {
            Self {
                min_hp_size: 3,
                temperature: 37.0,

                stack: *params.stack_en37,
                mismatch_hairpin: *params.mismatch_hairpin_en37,
                mismatch_interior: *params.mismatch_interior_en37,
                mismatch_interior_1n: *params.mismatch_interior_1n_en37,
                mismatch_interior_23: *params.mismatch_interior_23_en37,
                mismatch_multi: *params.mismatch_multi_en37,
                mismatch_exterior: *params.mismatch_exterior_en37,

                dangle5: *params.dangle5_en37,
                dangle3: *params.dangle3_en37,

                int11: *params.int11_en37,
                int21: *params.int21_en37,
                int22: *params.int22_en37,

                hairpin: *params.hairpin_en37,
                bulge: *params.bulge_en37,
                interior: *params.interior_en37,

                duplex_init: params.duplex_init_en37,
        
                closing_penalty: *params.closing_en37,

                lxc: params.lxc,

                ninio: params.ninio_en37,
                ninio_max: params.ninio_max,

                ml_base: params.ml_base_en37,
                ml_closing: params.ml_closing_en37,
                ml_intern: params.ml_intern_en37,

                triloops: params.triloops_en37.to_vec(),
                tetraloops: params.tetraloops_en37.to_vec(),
                hexaloops: params.hexaloops_en37.to_vec(),
            }
        } else {
            let kelvin = celsius + ZERO_CELSIUS_KELVIN;
            let scale = kelvin / (T_REF_CELSIUS + ZERO_CELSIUS_KELVIN);
            Self {
                min_hp_size: 3,
                temperature: celsius,

                stack: rescale_params!(stack, params, scale),
                mismatch_hairpin: rescale_params!(mismatch_hairpin, params, scale),
                mismatch_interior: rescale_params!(mismatch_interior, params, scale),
                mismatch_interior_1n: rescale_params!(mismatch_interior_1n, params, scale),
                mismatch_interior_23: rescale_params!(mismatch_interior_23, params, scale),
                mismatch_multi: rescale_params!(mismatch_multi, params, scale),
                mismatch_exterior: rescale_params!(mismatch_exterior, params, scale),

                dangle5: rescale_params!(dangle5, params, scale),
                dangle3: rescale_params!(dangle3, params, scale),

                int11: rescale_params!(int11, params, scale),
                int21: rescale_params!(int21, params, scale),
                int22: rescale_params!(int22, params, scale),

                hairpin: rescale_params!(hairpin, params, scale),
                bulge: rescale_params!(bulge, params, scale),
                interior: rescale_params!(interior, params, scale),

                duplex_init: rescale_param!(duplex_init, params, scale),

                closing_penalty: rescale_params!(closing, params, scale),

                lxc: params.lxc * scale,

                ninio: rescale_param!(ninio, params, scale),
                ninio_max: params.ninio_max,

                ml_base: rescale_param!(ml_base, params, scale),
                ml_closing: rescale_param!(ml_closing, params, scale),
                ml_intern: rescale_param!(ml_intern, params, scale),

                triloops: rescale_param!(triloops, params, scale),
                tetraloops: rescale_param!(tetraloops, params, scale),
                hexaloops: rescale_param!(hexaloops, params, scale),
            }
        }
    }

    /// The energy of each loop of a structure: `(0, description, energy)`, 
    /// in the order of the decomposition. They sum to the energy of the 
    /// structure. (Level 0, as in [`crate::Fuzzyfold::energy_terms`].)
    pub fn energy_terms<T: RingDecomposition>(&self, sequence: &[Base], structure: &T
    ) -> Result<Vec<(usize, String, i32)>, EnergyError> {
        let mut terms = Vec::new();
        self.sum_of_rings(sequence, structure, &mut |r, e| terms.push((0, r.to_string(), e)))?;
        Ok(terms)
    }

    /// The sum of the rings; each ring is also passed to `note`.
    fn sum_of_rings<T: RingDecomposition, N: FnMut(&Ring, i32)>(
        &self, sequence: &[Base], structure: &T, note: &mut N,
    ) -> Result<i32, EnergyError> {
        check_circular(sequence)?;
        let mut total = Ok(0);
        structure.for_each_ring(|r| {
            let Ok(sum) = total else { return };
            let en = match r {
                Ring::Loop(l) => self.energy_of_nested_loop(sequence, l),
                Ring::Pseudoknot(pk) => self.energy_of_pseudoknot(sequence, pk),
            };
            total = en.map(|en| {
                note(r, en);
                sum + en
            });
        });
        total
    }

    fn hairpin_bonus(&self, seq: &[Base]) -> Option<i32> {
        let table = match seq.len() {
            5 => &self.triloops,
            6 => &self.tetraloops,
            8 => &self.hexaloops,
            _ => return None,
        };
        table
            .iter()
            .find(|e| e.seq == seq)
            .map(|e| e.val)
    }

    /// The initiation of a hairpin loop with n unpaired bases: the table up
    /// to 30, extrapolated with lxc beyond.
    pub(super) fn hairpin_initiation(&self, n: usize) -> Result<i32, EnergyError> {
        if n < self.min_hp_size {
            return Err(EnergyError::HairpinTooSmall {
                size: n,
                min: self.min_hp_size,
            });
        }
        Ok(if n <= 30 {
            self.hairpin[n]
        } else {
            self.hairpin[30] + (self.lxc * ((n as f64) / 30.).ln()) as i32
        })
    }

    pub(super) fn hairpin(&self, seq: &[Base]) -> Result<i32, EnergyError> {
        let n = seq.len() - 2;
        let mut en = self.hairpin_initiation(n)?;

        // Special hairpin energies
        if let Some(en) = self.hairpin_bonus(seq) {
            return Ok(en);
        }

        let fb_closing = PairTypeRNA::from_fallback((seq[0], *seq.last().unwrap()));
        if !fb_closing.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }

        let closing_pair = PairTypeRNA::from((seq[0], *seq.last().unwrap()));

        if n == 3 {
            en += self.closing_penalty[closing_pair as usize];
        } else if n > 3 {
            en += self.mismatch_hairpin
                [fb_closing as usize]
                [seq[1].canonical_rna_index()]
                [seq[n].canonical_rna_index()]
                + self.closing_penalty[closing_pair as usize]
                - self.closing_penalty[fb_closing as usize];
        }

        Ok(en)
    }

    pub(super) fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> Result<i32, EnergyError> {
        let fb_outer = PairTypeRNA::from_fallback((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let fb_inner = PairTypeRNA::from_fallback((*revseq.first().unwrap(), *fwdseq.last().unwrap()));
        if !fb_outer.can_pair() || !fb_inner.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }

        let outer = PairTypeRNA::from((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));
        let pg1 = self.closing_penalty[outer as usize] -  self.closing_penalty[fb_outer as usize];
        let pg2 = self.closing_penalty[inner as usize] -  self.closing_penalty[fb_inner as usize];

        
        let res = match (fwdseq.len(), revseq.len()) {
            (2, 2) => {
                self.stack[outer as usize][inner as usize]
                    .ok_or(EnergyError::UnsupportedStacking { outer, inner })?
            },
            (3, 2) | (2, 3) => { //NOTE: SpecialC if C adjacent to paired C missing!
                self.bulge[1] + 
                    self.stack[outer as usize][inner as usize]
                    .ok_or(EnergyError::UnsupportedStacking { outer, inner })?
            },
            (3, 3) => 
                self.int11[fb_outer as usize][fb_inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [revseq[1].canonical_rna_index()] + pg1 + pg2,
            (3, 4) => 
                self.int21
                [fb_outer as usize][fb_inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [revseq[1].canonical_rna_index()]
                [revseq[2].canonical_rna_index()] + pg1 + pg2,
            (4, 3) => 
                self.int21
                [fb_inner as usize][fb_outer as usize]
                [revseq[1].canonical_rna_index()]
                [fwdseq[1].canonical_rna_index()]
                [fwdseq[2].canonical_rna_index()] + pg1 + pg2,
            (4, 4) => 
                self.int22
                [fb_outer as usize][fb_inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [fwdseq[2].canonical_rna_index()]
                [revseq[1].canonical_rna_index()]
                [revseq[2].canonical_rna_index()] + pg1 + pg2,
            (l, 2) | (2, l) => { // General Bulge case
                let n = l - 2;
                let ru_pg1 = self.closing_penalty[fb_outer as usize];
                let ru_pg2 = self.closing_penalty[fb_inner as usize];
                if n <= 30 {
                    self.bulge[n] + ru_pg1 + ru_pg2 + pg1 + pg2
                } else {
                    self.bulge[30] + ru_pg1 + ru_pg2 + pg1 + pg2
                    + (self.lxc * ((n as f64) / 30.).ln()) as i32
                }
            },
            (l, 3) | (3, l) => { // 1-n interior looop
                let mut en = 
                    self.mismatch_interior_1n
                    [fb_outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[revseq.len() - 2].canonical_rna_index()] +
                    self.mismatch_interior_1n
                    [fb_inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[fwdseq.len() - 2].canonical_rna_index()] + pg1 + pg2;

                en += self.ninio_max.min(
                    (l - 3) as i32 * self.ninio);

                let n = l - 1; 
                if n <= 30 {
                    en + self.interior[n]
                } else {
                    en + self.interior[30]
                       + (self.lxc * ((n as f64) / 30.).ln()) as i32
                }
            }
            (5, 4) | (4, 5) => { // 2-3 interior looop
                let mut en = 
                    self.mismatch_interior_23
                    [fb_outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[revseq.len() - 2].canonical_rna_index()] +
                    self.mismatch_interior_23
                    [fb_inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[fwdseq.len() - 2].canonical_rna_index()] + pg1 + pg2;
                en += self.ninio;
                en += self.interior[5];
                en
            }
            (lfwd, lrev) => { 
                let mut en = self.mismatch_interior
                    [fb_outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[lrev - 2].canonical_rna_index()] +
                    self.mismatch_interior
                    [fb_inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[lfwd - 2].canonical_rna_index()] + pg1 + pg2;

                let asy = (lfwd as isize - lrev as isize).abs() as i32;
                en += self.ninio_max.min(asy * self.ninio);
 
                let n = lfwd + lrev - 4; 
                if n <= 30 {
                    en + self.interior[n]
                } else {
                    en + self.interior[30]
                       + (self.lxc * ((n as f64) / 30.).ln()) as i32
                }
            }
        };
        Ok(res)
    }

    /// Multiloop energy from the segments between its helix ends, in loop
    /// order: the first segment starts and the last segment ends with the
    /// closing pair.
    pub(super) fn multibranch<S: AsRef<[Base]>>(&self, segments: impl IntoIterator<Item = S>
    ) -> Result<i32, EnergyError> {
        let mut segments = segments.into_iter();
        let first = segments.next().expect("A multiloop has segments.");
        let (mut en, mut n, mut m): (i32, usize, usize) = (0, 1, first.as_ref().len() - 2);
        let mut last: Option<S> = None;
        for s in segments {
            let a = last.as_ref().map_or(first.as_ref(), |l| l.as_ref());
            en += self.loop_helix_end(a, s.as_ref(), &self.mismatch_multi)?;
            (n, m) = (n + 1, m + s.as_ref().len() - 2);
            last = Some(s);
        }
        let last = last.as_ref().map_or(first.as_ref(), |l| l.as_ref());
        // For warning purposes only.
        if !PairTypeRNA::from((first.as_ref()[0], *last.last().unwrap())).can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }
        en += self.loop_helix_end(last, first.as_ref(), &self.mismatch_multi)?;
        Ok(en + self.ml_base * m as i32 + self.ml_closing + self.ml_intern * n as i32)
    }

    /// Exterior loop energy from the segments between its helix ends.
    pub(super) fn exterior<S: AsRef<[Base]>>(&self, segments: impl IntoIterator<Item = S>
    ) -> Result<i32, EnergyError> {
        let mut en = 0;
        let mut last: Option<S> = None;
        for s in segments {
            if let Some(a) = &last {
                en += self.loop_helix_end(a.as_ref(), s.as_ref(), &self.mismatch_exterior)?;
            }
            last = Some(s);
        }
        Ok(en)
    }

    /// The helix end between two consecutive segments of a loop: the last
    /// base of `a` and the first of `b`, flanked by their neighbors in the
    /// loop (none at a strand end).
    fn loop_helix_end(&self, a: &[Base], b: &[Base], mismatch: &MismatchParams) -> Result<i32, EnergyError> {
        let d5 = a.len().checked_sub(2).map(|d| a[d]);
        self.helix_end((*a.last().unwrap(), b[0]), d5, b.get(1).copied(), mismatch)
    }

    /// Terminal penalty plus mismatch or dangle energy of a helix end, as seen
    /// from the loop: `(p5, p3)` are the bases of the pair in loop order, `d5`
    /// and `d3` its neighbors in the loop.
    pub(super) fn helix_end(&self, 
        (p5, p3): (Base, Base), 
        d5: Option<Base>, 
        d3: Option<Base>, 
        mismatch: &MismatchParams,
    ) -> Result<i32, EnergyError> {
        let pair = PairTypeRNA::from((p5, p3));
        if !pair.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }
        Ok(self.closing_penalty[pair as usize] + self.dangles((p5, p3), d5, d3, mismatch))
    }

    /// Mismatch or dangle energy of a helix end, without the terminal penalty.
    pub(super) fn dangles(&self, 
        (p5, p3): (Base, Base), 
        d5: Option<Base>, 
        d3: Option<Base>, 
        mismatch: &MismatchParams,
    ) -> i32 {
        let fb_pair = PairTypeRNA::from_fallback((p5, p3));
        //NOTE: This does not take the minimum over all options, it always
        // prefers terminal mismatch over single dangling.
        // Also, in contrast to other mismatch conributions, it is also corrected for GU??
        match (d5, d3) { 
            (Some(b5), Some(b3)) => 
                mismatch[fb_pair as usize][b5.canonical_rna_index()][b3.canonical_rna_index()],
            (Some(b5), None) => 
                self.dangle5[fb_pair as usize][b5.canonical_rna_index()].min(0),
            (None, Some(b3)) => 
                self.dangle3[fb_pair as usize][b3.canonical_rna_index()].min(0),
            _ => 0,
        }
    }
}

/// Circular strands are evaluated alone: one strand, closed by the strand
/// wrap at its end. Anything else with a strand wrap is not supported in
/// ViennaRNA (as RNAeval -c).
fn check_circular(sequence: &[Base]) -> Result<(), EnergyError> {
    match sequence.iter().position(|&b| b == Base::SW) {
        None => Ok(()),
        Some(w) if w + 1 == sequence.len() && !sequence.contains(&Base::SB) => Ok(()),
        Some(_) => Err(EnergyError::CircularStrand),
    }
}

impl EnergyModel for ViennaRNA {
 
    fn temperature(&self) -> f64 {
        self.temperature
    }

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        PairTypeRNA::from((b1, b2)).can_pair()
    }

    fn min_hairpin_size(&self) -> usize { self.min_hp_size }

    fn energy_of_structure<T: RingDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> Result<i32, EnergyError> {
        self.sum_of_rings(sequence, structure, &mut |_, _| ())
    }

    fn energy_of_pseudoknot(&self, _sequence: &[Base], _pseudoknot: &Pseudoknot
    ) -> Result<i32, EnergyError> {
        Err(EnergyError::UnsupportedPseudoknot)
    }

    fn energy_of_nested_loop(&self, sequence: &[Base], nested: &NestedLoop
    ) -> Result<i32, EnergyError> {
        let en = self.nested_loop(sequence, nested)?;
        // Joins two strands: the 5' end after the 3' end.
        match nested.gap() {
            Gap::Open { five, three } if five > three => Ok(en + self.duplex_init),
            _ => Ok(en),
        }
    }
}

impl ViennaRNA {
    /// The energy of a loop, without the duplex initiation of a loop
    /// that joins two strands.
    pub(super) fn nested_loop(&self, sequence: &[Base], nested: &NestedLoop
    ) -> Result<i32, EnergyError> {
        let pairs = nested.pairs();
        match nested.gap() {
            // Stacks and hairpins without building segments.
            Gap::None => match pairs {
                [(j, i)] => self.hairpin(&sequence[*i..=*j]),
                [(p, q), (j, i)] => self.interior(&sequence[*i..=*p], &sequence[*q..=*j]),
                _ => self.multibranch(nested.segments(sequence)),
            },
            Gap::Open { .. } => self.exterior(nested.segments(sequence)),
            Gap::Wrap { .. } if pairs.is_empty() => Ok(0), // The open circle.
            Gap::Wrap { .. } => match nested.segments(sequence).collect::<Vec<_>>().as_slice() {
                [hairpin] => self.hairpin(hairpin),
                [left, right] => self.interior(left, right),
                segments => self.multibranch(segments),
            },
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::MultiPairTable;
    use crate::NucleotideVec;

    macro_rules! assert_hp {
        ($model:expr, $seq:expr, $val:expr) => {
            assert_eq!(
                $model
                .hairpin(&NucleotideVec::try_from_rna($seq).unwrap())
                .unwrap(),
                $val
            );
        };
    }

    macro_rules! assert_il {
        ($model:expr, $seq1:expr, $seq2:expr, $val:expr) => {
            assert_eq!(
                $model
                .interior(&NucleotideVec::try_from($seq1).unwrap(),
                          &NucleotideVec::try_from($seq2).unwrap())
                .unwrap(),
                $val
            );
        };
    }

    #[test]
    fn test_vrna_hairpin_evaluation() {
        let model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
        assert_hp!(model, "GAAAC", 540);
        assert_hp!(model, "CCGAGG", 350);
        assert_hp!(model, "CCAAGG", 330);
        assert_hp!(model, "CAAGG", 540);
        assert_hp!(model, "CAAAG", 540);
        assert_hp!(model, "AAAAU", 590);
        assert_hp!(model, "GAAAU", 590);
        assert_hp!(model, "CAAAAG", 410);
        assert_hp!(model, "ACCCU", 590);
        assert_hp!(model, "GCCCCC", 490);
        assert_hp!(model, "AAAAAU", 530);
        assert_hp!(model, "GAAAAU", 580);
        assert_hp!(model, "ACCCCU", 540);
        assert_hp!(model, "ACCCCCU", 550);
        assert_hp!(model, "AAAAAAU", 540);
        assert_hp!(model, "AAAAAAAU", 510);
        assert_hp!(model, "AAAAAAAAAAU", 610);
        assert_hp!(model, &format!("C{}G", "A".repeat(30)), 620);
        assert_hp!(model, &format!("C{}G", "A".repeat(31)), 623);
        assert_hp!(model, &format!("C{}G", "A".repeat(32)), 626);
        assert_hp!(model, &format!("C{}G", "A".repeat(33)), 630);
        assert_hp!(model, &format!("C{}G", "A".repeat(34)), 633);
        assert_hp!(model, &format!("C{}G", "A".repeat(35)), 636);
    }

    #[test]
    fn test_vrna_stacking_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CG", "CG", -240);
        assert_il!(model, "AC", "GU", -220);
        assert_il!(model, "GU", "AC", -220);
    }

    #[test]
    fn test_vrna_int11_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CCG", "CGG", 50);
        assert_il!(model, "CAG", "CAG", 90);
        assert_il!(model, "ACU", "AAU", 190);
        assert_il!(model, "GCU", "AUC", 120);
        assert_il!(model, "GCU", "AGC", 120);
    }

    #[test]
    fn test_vrna_int21_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CACG", "CGG", 110);
        assert_il!(model, "CAAG", "CAG", 230);
        assert_il!(model, "AACU", "AAU", 370);
        assert_il!(model, "GACU", "AUC", 300);
        assert_il!(model, "GACU", "AGC", 300);
        assert_il!(model, "CGG", "CACG", 110);
        assert_il!(model, "CAG", "CAAG", 230);
        assert_il!(model, "AAU", "AACU", 370);
        assert_il!(model, "AUC", "GACU", 300);
        assert_il!(model, "AGC", "GACU", 300);
    }

    #[test]
    fn test_vrna_bulge_1_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CAG", "CG", 140);
        assert_il!(model, "AAU", "AU", 270);
        assert_il!(model, "GAU", "AC", 160);
        assert_il!(model, "CCG", "CG", 140);
        assert_il!(model, "ACU", "AU", 270);
        assert_il!(model, "GCU", "AC", 160);
        assert_il!(model, "CG", "CAG", 140);
        assert_il!(model, "AU", "AAU", 270);
        assert_il!(model, "AC", "GAU", 160);
        assert_il!(model, "CG", "CCG", 140);
        assert_il!(model, "AU", "ACU", 270);
        assert_il!(model, "AC", "GCU", 160);
    }

    #[test]
    fn test_vrna_bulge_2_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CAAG", "CG", 280);
        assert_il!(model, "AAAU", "AU", 380);
        assert_il!(model, "GAAU", "AC", 330);
        assert_il!(model, "CCAG", "CG", 280);
        assert_il!(model, "ACAU", "AU", 380);
        assert_il!(model, "GCAU", "AC", 330);
        assert_il!(model, "CG", "CAAG", 280);
        assert_il!(model, "AU", "AAAU", 380);
        assert_il!(model, "AC", "GAAU", 330);
        assert_il!(model, "CG", "CCAG", 280);
        assert_il!(model, "AU", "ACAU", 380);
        assert_il!(model, "AC", "GCAU", 330);
    }

    #[test]
    fn test_vrna_interior_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "ACA", "UGAAU", 370);
        assert_il!(model, "ACAA", "UGAAU", 290);
        assert_il!(model, "GUAGU", "AGGC", 260);
        assert_il!(model, "AUAGU", "AGGU", 330);
        assert_il!(model, "GGC", "GUGC", 110);
    }

    #[test]
    fn test_vrna_bulge_n_evaluation() {
        let model = ViennaRNA::default();
        assert_il!(model, "CAAAAAAG", "CG", 440);
        assert_il!(model, "CAAAAAAAAG", "CG", 470);
    }

    #[test]
    fn test_vrna_multibranch() {
        let model = ViennaRNA::default();
        let seg1 = &NucleotideVec::try_from("GAAC").unwrap();
        let seg2 = &NucleotideVec::try_from("GAC").unwrap();
        let seg3 = &NucleotideVec::try_from("GAAAC").unwrap();
        let energy = model.multibranch([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, 330);
        let seg1 = &NucleotideVec::try_from("GAAC").unwrap();
        let seg2 = &NucleotideVec::try_from("GAC").unwrap();
        let seg3 = &NucleotideVec::try_from("GAAAAAAAAAAAAAAAAAAC").unwrap();
        let energy = model.multibranch([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, 330);
    }

    #[test]
    fn test_vrna_exterior_single_branch() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..]]).unwrap();
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::try_from("UG").unwrap();
        let seg2 = &NucleotideVec::try_from("CU").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..]]).unwrap();
        assert_eq!(energy, -120); 

        let seg1 = &NucleotideVec::try_from("G").unwrap();
        let seg2 = &NucleotideVec::try_from("CU").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..]]).unwrap();
        assert_eq!(energy, -120);
 
        let seg1 = &NucleotideVec::try_from("UG").unwrap();
        let seg2 = &NucleotideVec::try_from("C").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..]]).unwrap();
        assert_eq!(energy, 0); 
    }

    #[test]
    fn test_vrna_exterior_two_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUG").unwrap();
        let seg3 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUUG").unwrap();
        let seg3 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUUG").unwrap();
        let seg3 = &NucleotideVec::try_from("C").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CG").unwrap();
        let seg3 = &NucleotideVec::try_from("CU").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, -290);

        let seg1 = &NucleotideVec::try_from("ACA").unwrap();
        let seg2 = &NucleotideVec::try_from("UGG").unwrap();
        let seg3 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..]]).unwrap();
        assert_eq!(energy, -130);
    }

    #[test]
    fn test_vrna_exterior_three_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUG").unwrap();
        let seg3 = &NucleotideVec::try_from("CUG").unwrap();
        let seg4 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..], &seg4[..]]).unwrap();
        assert_eq!(energy, -360);

        let seg1 = &NucleotideVec::try_from("AUG").unwrap();
        let seg2 = &NucleotideVec::try_from("CUG").unwrap();
        let seg3 = &NucleotideVec::try_from("UUG").unwrap();
        let seg4 = &NucleotideVec::try_from("CUG").unwrap();
        let energy = model.exterior([&seg1[..], &seg2[..], &seg3[..], &seg4[..]]).unwrap();
        assert_eq!(energy, -240);
    }

    macro_rules! assert_eos {
        ($model:expr, $seq:expr, $dbr:expr, $val:expr) => {
            assert_eq!(
                $model
                .energy_of_structure(&NucleotideVec::try_from($seq).unwrap(),
                          &MultiPairTable::try_from($dbr).unwrap()).unwrap(),
                $val
            );
        };
    }
 
    #[test]
    fn test_evaluations() {
        let model = ViennaRNA::default();

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eos!(model, seq, dbr, e37);

        let seq = "ACGUUAAAGACGU";
        let dbr = "(((((...)))))";
        let e37 = -170;
        assert_eos!(model, seq, dbr, e37);

        let seq = "AGACGACAAGGUUGAAUCGC";
        let dbr = ".(.(((.(....)...))))";
        let e37 = 420;
        assert_eos!(model, seq, dbr, e37);

        let seq = "GAGUAGUGGAACCAGGCUAU";
        let dbr = ".((...((....))..))..";
        let e37 = 190;
        assert_eos!(model, seq, dbr, e37);

        let seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACC";
        let dbr = "...........(.(((((........)))))..)......";
        let e37 = -210;
        assert_eos!(model, seq, dbr, e37);
    }
 
    macro_rules! assert_meos {
        ($model:expr, $seq:expr, $dbr:expr, $val:expr) => {
            assert_eq!(
                $model
                .energy_of_structure(&NucleotideVec::try_from($seq).unwrap(),
                          &MultiPairTable::try_from($dbr).unwrap()).unwrap(),
                $val
            );
        };
    }

    #[test]
    fn test_multi_evaluations() {
        let model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_meos!(model, seq, dbr, e37);

        let seq = "GA+AAAC";
        let dbr = "(.+...)";
        let e37 = 300;
        assert_meos!(model, seq, dbr, e37);
 
        let seq = "AAAC+GA";
        let dbr = "...(+).";
        let e37 = 300;
        assert_meos!(model, seq, dbr, e37);

        let seq = "GAA+AAC";
        let dbr = "(..+..)";
        let e37 = 300;
        assert_meos!(model, seq, dbr, e37);

        let seq = "GC+UUUUAGU+AU+AC";
        let dbr = "((+(...)).+..+.)";
        let e37 = 1140;
        assert_meos!(model, seq, dbr, e37);

        let seq = "GC&UUUUAGU&AGAAACU&AGAAACU&AC";
        let dbr = "((&(...)).&.(...).&.(...).&.)";
        let e37 = 2020;
        assert_meos!(model, seq, dbr, e37);
    }

    #[test]
    fn parses_inosine() {
        assert_eq!(Base::try_from('I').unwrap(), Base::IN);
    }

    #[test]
    fn inosine_in_rna() {
        assert!(NucleotideVec::try_from_rna("AGI").is_ok());
    }

     #[test]
    fn inosine_canonical_rna_index() {
        assert_eq!(Base::IN.canonical_rna_index(), Base::G as usize);
    }
    
    #[test]
    fn test_circular_evaluations() {
        // Energies of ViennaRNA 2.6.4 for circular RNA (RNAeval -c). The
        // exterior loop of the circle is the loop through the strand wrap.
        let model = ViennaRNA::default();
        // The open circle.
        assert_eos!(model, "ACGUACGU^", "........^", 0);
        // Hairpins through the wrap, the first with the special hairpin CUUCGG.
        assert_eos!(model, "CGGGACAAAAGUCCUU^", "..((((....))))..^", -10);
        assert_eos!(model, "ACGCUGUCUGAGACUAGAAGACAGAUAGU^", "...((((((.........)))))).....^", -180);
        // A stack and an interior loop through the wrap.
        assert_eos!(model, "UGUCUGAGACUAGAAGACAGAUAGUACGC^", "(((((.........)))))(........)^", -180);
        assert_eos!(model, "AGACUAGAAGACAGAUAGUACGCUGUCUG^", "...((...)).(((((((....)))))))^", -50);
        // A multiloop through the wrap.
        assert_eos!(model, "AAUGCUAAGACAUUUCCCUUCAGGGGGGGCUCCCCCGC^",
                           "((((......)))).(((....)))((((...))))..^", -290);
    }

    #[test]
    fn test_circular_strands_not_supported() {
        let model = ViennaRNA::default();
        let energy = |seq: &str, db: &str| model.energy_of_structure(
            &NucleotideVec::try_from(seq).unwrap(), &MultiPairTable::try_from(db).unwrap());
        // The loop through the wrap is a hairpin without unpaired bases.
        assert!(matches!(energy("GGGAAACCC^", "(((...)))^"), Err(EnergyError::HairpinTooSmall { .. })));
        // A circular strand together with another strand.
        assert!(matches!(energy("GGGAAACCC^AAA", "(((...)))^..."), Err(EnergyError::CircularStrand)));
    }

    #[test]
    fn rejects_pseudoknots() {
        let vrna = ViennaRNA::default();
        let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let pk = MultiPairTable::try_from("[[[....(((((]]]......)))))").unwrap();
        assert!(matches!(vrna.energy_of_structure(&seq, &pk), Err(EnergyError::UnsupportedPseudoknot)));
        assert!(matches!(vrna.energy_terms(&seq, &pk), Err(EnergyError::UnsupportedPseudoknot)));

        // Without pseudoknots, both models agree on the same parameters.
        let nested = MultiPairTable::try_from(".......(((((.........)))))").unwrap();
        let ff = crate::Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
        assert_eq!(vrna.energy_of_structure(&seq, &nested).unwrap(), 
            ff.energy_of_structure(&seq, &nested).unwrap());
    }
}
