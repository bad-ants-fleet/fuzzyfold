use log::info;
use colored::*; 
use paste::paste;

use crate::Base;
use crate::PairTypeRNA;
use crate::EnergyModel;
use crate::EnergyError;
use crate::parameters::*;
use crate::LoopDecomposition;
use crate::NearestNeighborLoop;
use crate::K0;
use crate::nucleotides::BCOUNT;

pub struct ViennaRNA {
    min_hp_size: usize,
    temperature: f64, //TODO: make this optional for fitted params?

    stack: ExtendedStackParams,
    mismatch_hairpin: MismatchParams,
    mismatch_interior: MismatchParams,
    mismatch_interior_1n: MismatchParams,
    mismatch_interior_23: MismatchParams,
    mismatch_multi: MismatchParams,
    mismatch_exterior: MismatchParams,
                                   
    dangle5: DangleParams,
    dangle3: DangleParams,
                                                                              
    int11: Int11Params,
    int21: Int21Params,
    int22: Int22Params,

    hairpin: LoopParams,
    bulge: LoopParams,
    interior: LoopParams,

    duplex_init: i32,
    terminal_ru: i32,
    lxc: f64,

    ninio: i32,
    ninio_max: i32,

    ml_base: i32,
    ml_closing: i32,
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
    pub fn from_fitted_params(params: &'static FittedParams) -> Self {
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
            terminal_ru: params.terminal_ru,
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
    pub fn from_thermo_params(params: &'static ThermoParams, celsius: f64) -> Self {
        if (celsius - T_REF).abs() < 1e-6 {
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
                terminal_ru: params.terminal_ru_en37,
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
            let kelvin = celsius + K0;
            let scale = kelvin / (T_REF + K0);
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
                terminal_ru: rescale_param!(terminal_ru, params, scale),
                lxc: params.lxc * celsius,

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

    fn hairpin(&self, seq: &[Base]) -> Result<i32, EnergyError> {
        let n = seq.len() - 2;

        if n < self.min_hp_size {
            return Err(EnergyError::HairpinTooSmall {
                size: n,
                min: self.min_hp_size,
            });
        }

        // Special hairpin energies
        if let Some(en) = self.hairpin_bonus(seq) {
            return Ok(en);
        }

        let closing = PairTypeRNA::from_fallback((seq[0], *seq.last().unwrap()));
        if !closing.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }

        // Initiation term
        let mut en = if n <= 30 {
            self.hairpin[n]
        } else {
            self.hairpin[30] + (self.lxc * ((n as f64) / 30.).ln()) as i32
        };

        if n == 3 && closing.is_ru() {
            en += self.terminal_ru;
        } else if n > 3 {
            en += self.mismatch_hairpin
                [closing as usize]
                [seq[1].canonical_rna_index()]
                [seq[n].canonical_rna_index()];
        }

        Ok(en)
    }

    fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> Result<i32, EnergyError> {
        let outer = PairTypeRNA::from_fallback((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from_fallback((*revseq.first().unwrap(), *fwdseq.last().unwrap()));
        if !outer.can_pair() || !inner.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }

        let res = match (fwdseq.len(), revseq.len()) {
            (2, 2) => {
                let outer = PairTypeRNA::from((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
                let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));
                self.stack[outer as usize][inner as usize]
            },
            (3, 2) | (2, 3) => { //NOTE: SpecialC if C adjacent to paired C missing!
                self.bulge[1] + 
                self.stack[outer as usize][inner as usize]},
            (3, 3) => 
                self.int11[outer as usize][inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [revseq[1].canonical_rna_index()],
            (3, 4) => 
                self.int21
                [outer as usize][inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [revseq[1].canonical_rna_index()]
                [revseq[2].canonical_rna_index()],
            (4, 3) => 
                self.int21
                [inner as usize][outer as usize]
                [revseq[1].canonical_rna_index()]
                [fwdseq[1].canonical_rna_index()]
                [fwdseq[2].canonical_rna_index()],
            (4, 4) => 
                self.int22
                [outer as usize][inner as usize]
                [fwdseq[1].canonical_rna_index()]
                [fwdseq[2].canonical_rna_index()]
                [revseq[1].canonical_rna_index()]
                [revseq[2].canonical_rna_index()],
            (l, 2) | (2, l) => { // General Bulge case
                let n = l - 2;
                let pg1 = if outer.is_ru() { self.terminal_ru } else { 0 };
                let pg2 = if inner.is_ru() { self.terminal_ru } else { 0 };
                if n <= 30 {
                    self.bulge[n] + pg1 + pg2
                } else {
                    self.bulge[30] + pg1 + pg2
                    + (self.lxc * ((n as f64) / 30.).ln()) as i32
                }
            },
            (l, 3) | (3, l) => { // 1-n interior looop
                let mut en = 
                    self.mismatch_interior_1n
                    [outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[revseq.len() - 2].canonical_rna_index()] +
                    self.mismatch_interior_1n
                    [inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[fwdseq.len() - 2].canonical_rna_index()];

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
                    [outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[revseq.len() - 2].canonical_rna_index()] +
                    self.mismatch_interior_23
                    [inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[fwdseq.len() - 2].canonical_rna_index()];
                en += self.ninio;
                en += self.interior[5];
                en
            }
            (lfwd, lrev) => { 
                let mut en = self.mismatch_interior
                    [outer as usize]
                    [fwdseq[1].canonical_rna_index()]
                    [revseq[lrev - 2].canonical_rna_index()] +
                    self.mismatch_interior
                    [inner as usize]
                    [revseq[1].canonical_rna_index()]
                    [fwdseq[lfwd - 2].canonical_rna_index()];

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

    fn multibranch(&self, segments: &[&[Base]]) -> Result<i32, EnergyError> {
        // For warning purposes only.
        let closing = PairTypeRNA::from((segments[0][0], *segments.last().unwrap().last().unwrap()));
        if !closing.can_pair() {
            return Err(EnergyError::InvalidClosingPair);
        }

        // Number of stems in the multiloop.
        let n = segments.len(); 

        let mut en = 0;
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = PairTypeRNA::from_fallback((*segments[i].last().unwrap(), segments[j][0]));
            if !pair.can_pair() {
                return Err(EnergyError::InvalidClosingPair);
            }
            if pair.is_ru() { 
                en += self.terminal_ru;
            }
            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j).and_then(|seg| seg.get(1));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let den = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    self.mismatch_multi
                    [pair as usize]
                    [b5.canonical_rna_index()]
                    [b3.canonical_rna_index()],
                (Some(&b5), None) => 
                    self.dangle5
                     [pair as usize]
                     [b5.canonical_rna_index()].min(0),
                (None, Some(&b3)) => 
                    self.dangle3
                    [pair as usize]
                    [b3.canonical_rna_index()].min(0),
                _ => 0,
            };
            en += den;
        }
 
        // Number of unpaired bases in the multiloop.
        let m: usize = segments.iter().map(|s| s.len() - 2).sum();
        Ok(en + self.ml_base * m as i32
           + self.ml_closing
           + self.ml_intern * n as i32)
    }

    fn exterior(&self, segments: &[&[Base]]) -> Result<i32, EnergyError> {
        let mut en = 0;
        let n = segments.len() - 1; 
        for i in 0..n {
            let j = i + 1; 
            
            let pair = PairTypeRNA::from_fallback((*segments[i].last().unwrap(), segments[j][0]));
            if !pair.can_pair() {
                return Err(EnergyError::InvalidClosingPair);
            }
            if pair.is_ru() { 
                en += self.terminal_ru;
            }

            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j).and_then(|seg| seg.get(1));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let den = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    self.mismatch_exterior
                    [pair as usize][b5.canonical_rna_index()][b3.canonical_rna_index()],
                (Some(&b5), None) => 
                    self.dangle5
                    [pair as usize][b5.canonical_rna_index()].min(0),
                (None, Some(&b3)) => 
                     self.dangle3
                    [pair as usize][b3.canonical_rna_index()].min(0),
                _ => 0,
            };
            en += den;
        }
        Ok(en)
    }
}

const CAN_PAIR: [[bool; BCOUNT]; BCOUNT] = {
    use Base::*;
    let mut table = [[false; BCOUNT]; BCOUNT];
    table[A as usize][U as usize] = true;
    table[U as usize][A as usize] = true;
    table[C as usize][G as usize] = true;
    table[G as usize][C as usize] = true;
    table[G as usize][U as usize] = true;
    table[U as usize][G as usize] = true;
    table[A as usize][PU as usize] = true;
    table[PU as usize][A as usize] = true;
    table[G as usize][PU as usize] = true;
    table[PU as usize][G as usize] = true;
    table
};

impl EnergyModel for ViennaRNA {
 
    fn temperature(&self) -> f64 {
        self.temperature
    }

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        CAN_PAIR[b1 as usize][b2 as usize]
    }

    fn min_hairpin_size(&self) -> usize { self.min_hp_size }

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> Result<i32, EnergyError> {
        let mut total = 0;
        structure.for_each_loop(|l| {
            let en = self.energy_of_loop(sequence, l).unwrap_or_else(|e| {
                panic!("Energy evaluation error: {:?} {:?} {:?}.", sequence, l, e);
            });
            total += en;
            info!("{:<41} {}", format!("{}:", l), format!("{:>6.2}", en as f64 / 100.).green());
        });
        Ok(total)
    }

    fn energy_of_loop(&self, sequence: &[Base], nn_loop: &NearestNeighborLoop
    ) -> Result<i32, EnergyError> {
        match nn_loop {
            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                self.hairpin(&sequence[*i as usize..=*j as usize])
            }
            NearestNeighborLoop::Interior { closing: (i, j), inner: (k, l) } => {
                let left = &sequence[*i as usize..=*k as usize];
                let right = &sequence[*l as usize..=*j as usize];
                self.interior(left, right)
            }
            NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                let mut slices: Vec<&[Base]> = Vec::with_capacity(branches.len() + 1);
                let mut start = *i as usize;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k as usize]);
                    start = l as usize;
                }
                slices.push(&sequence[start..=*j as usize]);
                self.multibranch(&slices)
            }
            NearestNeighborLoop::Exterior { ends: (p5, p3), branches  } => {
                let mut slices: Vec<&[Base]> = Vec::with_capacity(branches.len() + 1);
                let mut p5 = *p5 as usize;
                for &(k, l) in branches {
                    slices.push(&sequence[p5..=k as usize]);
                    p5 = l as usize;
                }
                slices.push(&sequence[p5..=(*p3 as usize)]);
                self.exterior(&slices)
            }
            NearestNeighborLoop::JointExterior { ends: (p5, p3), branches  } => {
                let mut slices: Vec<&[Base]> = Vec::with_capacity(branches.len() + 1);
                let mut branches = branches.clone();

                debug_assert!(!branches.is_empty());
                branches.rotate_left(1);
                let last = branches.len() - 1;
                let (i, j) = branches[last];
                branches[last] = (j, i);
                while let Some(&(i, _)) = branches.first() {
                    if i > *p3 { break; }
                    branches.rotate_left(1);
                }

                let mut p5 = *p5 as usize;
                for (k, l) in branches {
                    slices.push(&sequence[p5..=k as usize]);
                    p5 = l as usize;
                }
                slices.push(&sequence[p5..=(*p3 as usize)]);
                self.exterior(&slices).map(|e| e + self.duplex_init)
            }
            NearestNeighborLoop::Disconnected { .. } => unreachable!("Must not evaluate disconnected loops.")
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_structure::MultiPairTable;
    use crate::NucleotideVec;

    macro_rules! assert_hp {
        ($model:expr, $seq:expr, $val:expr) => {
            assert_eq!(
                $model
                .hairpin(&NucleotideVec::from_lossy($seq))
                .unwrap(),
                $val
            );
        };
    }

    macro_rules! assert_il {
        ($model:expr, $seq1:expr, $seq2:expr, $val:expr) => {
            assert_eq!(
                $model
                .interior(&NucleotideVec::from_lossy($seq1),
                          &NucleotideVec::from_lossy($seq2))
                .unwrap(),
                $val
            );
        };
    }

    #[test]
    fn test_vrna_hairpin_evaluation() {
        let model = ViennaRNA::default();
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
        let seg1 = &NucleotideVec::from_lossy("GAAC");
        let seg2 = &NucleotideVec::from_lossy("GAC");
        let seg3 = &NucleotideVec::from_lossy("GAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, 330);
        let seg1 = &NucleotideVec::from_lossy("GAAC");
        let seg2 = &NucleotideVec::from_lossy("GAC");
        let seg3 = &NucleotideVec::from_lossy("GAAAAAAAAAAAAAAAAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, 330);
    }

    #[test]
    fn test_vrna_exterior_single_branch() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2]).unwrap();
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::from_lossy("UG");
        let seg2 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2]).unwrap();
        assert_eq!(energy, -120); 

        let seg1 = &NucleotideVec::from_lossy("G");
        let seg2 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2]).unwrap();
        assert_eq!(energy, -120);
 
        let seg1 = &NucleotideVec::from_lossy("UG");
        let seg2 = &NucleotideVec::from_lossy("C");
        let energy = model.exterior(&[seg1, seg2]).unwrap();
        assert_eq!(energy, 0); 
    }

    #[test]
    fn test_vrna_exterior_two_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUUG");
        let seg3 = &NucleotideVec::from_lossy("C");
        let energy = model.exterior(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CG");
        let seg3 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, -290);

        let seg1 = &NucleotideVec::from_lossy("ACA");
        let seg2 = &NucleotideVec::from_lossy("UGG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]).unwrap();
        assert_eq!(energy, -130);
    }

    #[test]
    fn test_vrna_exterior_three_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let seg4 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]).unwrap();
        assert_eq!(energy, -360);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("UUG");
        let seg4 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]).unwrap();
        assert_eq!(energy, -240);
    }

 
    #[test]
    fn test_evaluations() {
        let model = ViennaRNA::default();

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "ACGUUAAAGACGU";
        let dbr = "(((((...)))))";
        let e37 = -170;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "AGACGACAAGGUUGAAUCGC";
        let dbr = ".(.(((.(....)...))))";
        let e37 = 420;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "GAGUAGUGGAACCAGGCUAU";
        let dbr = ".((...((....))..))..";
        let e37 = 190;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACC";
        let dbr = "...........(.(((((........)))))..)......";
        let e37 = -210;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")).unwrap(), e37);
    }
 
    #[test]
    fn test_multi_evaluations() {
        let model = ViennaRNA::default();

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let fseq = "GA+AAAC";
        let fdbr = "(.+...)";
 
        let rseq = "AAAC+GA";
        let rdbr = "...(+).";
        assert_eq!(
            model.energy_of_structure(
                &NucleotideVec::from_lossy(fseq), 
                &MultiPairTable::try_from(fdbr).expect("valid")).unwrap(), 
            model.energy_of_structure(
                &NucleotideVec::from_lossy(rseq), 
                &MultiPairTable::try_from(rdbr).expect("valid")).unwrap()
        );

        let seq = "GAA+AAC";
        let dbr = "(..+..)";
        let e37 = 300;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "GC+UUUUAGU+AU+AC";
        let dbr = "((+(...)).+..+.)";
        let e37 = 1140;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")).unwrap(), e37);

        let seq = "GC&UUUUAGU&AGAAACU&AGAAACU&AC";
        let dbr = "((&(...)).&.(...).&.(...).&.)";
        let e37 = 2020;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")).unwrap(), e37);
 
    }

}

