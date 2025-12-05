use log::info;
use colored::*; 

use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

use crate::parameters::TURNER_2004;
use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;
use crate::PairTypeRNA;
use crate::nn_parsing::EnergyTables;
use crate::nn_parsing::ParamError;
use crate::EnergyModel;
use crate::K0;

/// The default ViennaRNA-v2.6 energy model.
///
/// The current implementation allows different
/// parameter files, but only the standard
/// settings (dangle=2, tetraloops, etc.)
///
/// Only single-stranded folding is supported.
///
pub struct ViennaRNA {
    min_hp_size: usize,
    temperature: f64,
    energy_tables: EnergyTables,
}

impl Default for ViennaRNA {
    fn default() -> Self {
        ViennaRNA::from_parameter_str(TURNER_2004)
            .expect("Built-in Turner 2004 parameter file must be valid")
    }
}

impl ViennaRNA {

    pub fn from_parameter_file<P: AsRef<Path>>(path: P) -> Result<Self, ParamError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Self::from_reader(reader)
    }

    pub fn from_parameter_str(s: &str) -> Result<Self, ParamError> {
        let cursor = std::io::Cursor::new(s);
        Self::from_reader(cursor)
    }

    fn from_reader<R: BufRead>(reader: R) -> Result<Self, ParamError> {
        let energy_tables = EnergyTables::from_reader(reader)?;
        Ok(ViennaRNA {
            min_hp_size: 3,
            temperature: 37.0,
            energy_tables,
        })
    }

    pub fn set_temperature(&mut self, temperature: f64) {
        if (self.temperature - temperature).abs() < f64::EPSILON {
            return;
        }

        let old_temp = self.temperature + K0;
        let new_temp = temperature + K0;
        let temp_change = new_temp / old_temp;
        self.temperature = temperature;
        self.energy_tables.rescale(temp_change);
    }

    fn hairpin(&self, seq: &[Base]) -> i32 {
        let et = &self.energy_tables;
        let n = seq.len() - 2;
        if n < self.min_hp_size {
            panic!("Invalid hairpin size {n}");
        }
        let closing = PairTypeRNA::new((seq[0], *seq.last().unwrap()));

        // Special hairpin energies
        if seq.len() <= 6 { 
            if let Some((en, _)) = et.hairpin_sequences.get(seq).copied() {
                return en;
            }
        }

        // Initiation terms
        let mut en = if n <= 30 {
            et.hairpin[n].expect("from file")
        } else {
            et.hairpin[30].expect("from file")
            + (self.energy_tables.misc.lxc * ((n as f64) / 30.).ln()) as i32
        };

        // NOTE: double check NN desrciption if this is indeed the correct way.
        if n == 3 && matches!(closing
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA
            | PairTypeRNA::NN) { 
            en += et.misc.terminal_ru_en37;
        } else if n > 3 {
            en += et.mismatch_hairpin
                [closing as usize]
                [seq[1] as usize]
                [seq[n] as usize].expect("from file");
        }
        en
    }

    fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> i32 {
        let outer = PairTypeRNA::new((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));

        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA | PairTypeRNA::NN);

        match (fwdseq.len(), revseq.len()) {
            (2, 2) => 
                self.energy_tables.stack[outer as usize][inner as usize]
                .expect("from file"),
            (3, 2) | (2, 3) => //NOTE: SpecialC if C adjacent to paired C missing!
                self.energy_tables.bulge[1].unwrap() + 
                self.energy_tables.stack[outer as usize][inner as usize]
                .expect("from file"),
            (3, 3) => 
                self.energy_tables.int11[outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize]
                .expect("from file"),
            (3, 4) => 
                self.energy_tables.int21
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize]
                [revseq[2] as usize]
                .expect("from file"),
            (4, 3) => 
                self.energy_tables.int21
                [inner as usize][outer as usize]
                [revseq[1] as usize][fwdseq[1] as usize]
                [fwdseq[2] as usize].expect("from file"),
            (4, 4) => 
                self.energy_tables.int22
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize][revseq[2] as usize]
                .expect("from file"),
            (l, 2) | (2, l) => { // General Bulge case
                let n = l - 2;
                let pg1 = if !is_ru_end(outer) { 0 } else {
                    self.energy_tables.misc.terminal_ru_en37
                };
                let pg2 = if !is_ru_end(inner) { 0 } else {
                    self.energy_tables.misc.terminal_ru_en37
                };
                if n <= 30 {
                    self.energy_tables.bulge[n].unwrap() + pg1 + pg2
                } else {
                    self.energy_tables.bulge[30].unwrap() + pg1 + pg2
                    + (self.energy_tables.misc.lxc * ((n as f64) / 30.).ln()) as i32
                }
            },
            (l, 3) | (3, l) => { // 1-n interior looop
                let mut en = 
                    self.energy_tables.mismatch_interior_1n
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[revseq.len() - 2] as usize].unwrap() +
                    self.energy_tables.mismatch_interior_1n
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[fwdseq.len() - 2] as usize].unwrap();

                en += self.energy_tables.ninio.max.min(
                    (l - 3) as i32 * self.energy_tables.ninio.en37);

                let n = l - 1; 
                if n <= 30 {
                    en + self.energy_tables.interior[n].unwrap()
                } else {
                    en + self.energy_tables.interior[30].unwrap() 
                       + (self.energy_tables.misc.lxc * ((n as f64) / 30.).ln()) as i32
                }
            }
            (5, 4) | (4, 5) => { // 2-3 interior looop
                let mut en = 
                    self.energy_tables.mismatch_interior_23
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[revseq.len() - 2] as usize].unwrap() +
                    self.energy_tables.mismatch_interior_23
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[fwdseq.len() - 2] as usize].unwrap();
                en += self.energy_tables.ninio.en37;
                en += self.energy_tables.interior[5].unwrap();
                en
            }
            (lfwd, lrev) => { 
                let mut en = self.energy_tables.mismatch_interior
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[lrev - 2] as usize].unwrap() +
                    self.energy_tables.mismatch_interior
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[lfwd - 2] as usize].unwrap();

                let asy = (lfwd as isize - lrev as isize).abs() as i32;
                en += self.energy_tables.ninio.max.min(
                    asy * self.energy_tables.ninio.en37);
 
                let n = lfwd + lrev - 4; 
                if n <= 30 {
                    en + self.energy_tables.interior[n].unwrap()
                } else {
                    en + self.energy_tables.interior[30].unwrap() 
                       + (self.energy_tables.misc.lxc * ((n as f64) / 30.).ln()) as i32
                }
            }
        }
    }

    fn multibranch(&self, segments: &[&[Base]]) -> i32 {
        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA 
            | PairTypeRNA::NN);

        // For warning purposes only.
        let _ = PairTypeRNA::new((segments[0][0], *segments.last().unwrap().last().unwrap()));

        let n = segments.len(); 

        let mut en = 0;
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            if is_ru_end(pair) { 
                en += self.energy_tables.misc.terminal_ru_en37;
            }
            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j).and_then(|seg| seg.get(1));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let den = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    self.energy_tables.mismatch_exterior
                    [pair as usize][b5 as usize][b3 as usize].unwrap(),
                (Some(&b5), None) => 
                    self.energy_tables.dangle5
                     [pair as usize][b5 as usize].unwrap().min(0),
                (None, Some(&b3)) => 
                    self.energy_tables.dangle3
                    [pair as usize][b3 as usize].unwrap().min(0),
                _ => 0,
            };
            en += den;
        }
 
        en + self.energy_tables.ml_params.base_en37 
           + self.energy_tables.ml_params.closing_en37
           + self.energy_tables.ml_params.intern_en37 * n as i32
    }

    fn exterior(&self, segments: &[&[Base]]) -> i32 {
        let is_ru_end = |pt| matches!(pt
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA | PairTypeRNA::NN);
 
        let mut en = 0;
        let n = segments.len() - 1; 
        for i in 0..n {
            let j = i + 1; 
            
            let pair = (*segments[i].last().unwrap(), segments[j][0]).into();
            if is_ru_end(pair) { 
                en += self.energy_tables.misc.terminal_ru_en37;
            }

            let d5 = segments.get(i)
                .and_then(|seg| seg.len().checked_sub(2).and_then(|d| seg.get(d)));
            let d3 = segments.get(j).and_then(|seg| seg.get(1));

            //NOTE: This does not take the minimum over all options, it always
            // prefers terminal mismatch over single dangling.
            let den = match (d5, d3) { 
                (Some(&b5), Some(&b3)) => 
                    self.energy_tables.mismatch_exterior
                    [pair as usize][b5 as usize][b3 as usize].unwrap(),
                (Some(&b5), None) => 
                    self.energy_tables.dangle5
                    [pair as usize][b5 as usize].unwrap().min(0),
                (None, Some(&b3)) => 
                     self.energy_tables.dangle3
                    [pair as usize][b3 as usize].unwrap().min(0),
                _ => 0,
            };
            en += den;
        }
        en
    }

}

impl EnergyModel for ViennaRNA {
 
    fn temperature(&self) -> f64 {
        self.temperature
    }

    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        matches!((b1, b2),
        (Base::A, Base::U) | (Base::U, Base::A) |
        (Base::G, Base::C) | (Base::C, Base::G) |
        (Base::G, Base::U) | (Base::U, Base::G))
    }

    fn min_hairpin_size(&self) -> usize { self.min_hp_size }

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> i32  {
        let mut total = 0;
        structure.for_each_loop(|l| {
            let en = self.energy_of_loop(sequence, l);
            total += en;
            info!("{:<41} {}", format!("{}:", l), format!("{:>6.2}", en as f64 / 100.).green());
        });
        total
    }

    fn energy_of_loop(&self, sequence: &[Base], nn_loop: &NearestNeighborLoop) -> i32 {

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
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = *i as usize;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k as usize]);
                    start = l as usize;
                }
                slices.push(&sequence[start..=*j as usize]);
                self.multibranch(&slices)
            }
            NearestNeighborLoop::Exterior { branches } => {
                let mut slices: Vec<&[Base]> = Vec::new();
                let mut start = 0;
                for &(k, l) in branches {
                    slices.push(&sequence[start..=k as usize]);
                    start = l as usize;
                }
                slices.push(&sequence[start..]);
                self.exterior(&slices)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use crate::NucleotideVec;

    #[test]
    fn test_vrna_hairpin_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAC")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CCGAGG")), 350);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CCAAGG")), 330);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAGG")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAAG")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAU")), 590);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAU")), 590);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAAAG")), 410);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("ACCCU")), 590);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GCCCCC")), 490);

        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAG")), 590);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAAC")), 590);

        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAU")), 530);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAAU")), 580);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("ACCCCU")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("ACCCCCU")), 550);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAU")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAU")), 510);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAAAAU")), 610);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAAAAA")), 660);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(30)))), 620);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(31)))), 623);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(32)))), 626);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(33)))), 630);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(34)))), 633);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(35)))), 636);
    }

    #[test]
    fn test_vrna_stacking_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CG")), -240);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GU")), -220);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GU"), &NucleotideVec::from_lossy("AC")), -220);
    }

    #[test]
    fn test_vrna_int11_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CCG"), &NucleotideVec::from_lossy("CGG")), 50);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAG"), &NucleotideVec::from_lossy("CAG")), 90);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACU"), &NucleotideVec::from_lossy("AAU")), 190);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCU"), &NucleotideVec::from_lossy("AUC")), 120);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCU"), &NucleotideVec::from_lossy("AGC")), 120);
    }

    #[test]
    fn test_vrna_int21_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CACG"), &NucleotideVec::from_lossy("CGG")), 110);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAG"), &NucleotideVec::from_lossy("CAG")), 230);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AACU"), &NucleotideVec::from_lossy("AAU")), 370);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GACU"), &NucleotideVec::from_lossy("AUC")), 300);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GACU"), &NucleotideVec::from_lossy("AGC")), 300);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CGG"), &NucleotideVec::from_lossy("CACG")), 110);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAG"), &NucleotideVec::from_lossy("CAAG")), 230);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AAU"), &NucleotideVec::from_lossy("AACU")), 370);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AUC"), &NucleotideVec::from_lossy("GACU")), 300);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AGC"), &NucleotideVec::from_lossy("GACU")), 300);
    }

    #[test]
    fn test_vrna_bulge_1_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAG"), &NucleotideVec::from_lossy("CG")), 140);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AAU"), &NucleotideVec::from_lossy("AU")), 270);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GAU"), &NucleotideVec::from_lossy("AC")), 160);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CCG"), &NucleotideVec::from_lossy("CG")), 140);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACU"), &NucleotideVec::from_lossy("AU")), 270);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCU"), &NucleotideVec::from_lossy("AC")), 160);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CAG")), 140);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AU"), &NucleotideVec::from_lossy("AAU")), 270);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GAU")), 160);
                                                                
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CCG")), 140);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AU"), &NucleotideVec::from_lossy("ACU")), 270);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GCU")), 160);
    }

    #[test]
    fn test_vrna_bulge_2_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAG"), &NucleotideVec::from_lossy("CG")), 280);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AAAU"), &NucleotideVec::from_lossy("AU")), 380);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GAAU"), &NucleotideVec::from_lossy("AC")), 330);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CCAG"), &NucleotideVec::from_lossy("CG")), 280);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACAU"), &NucleotideVec::from_lossy("AU")), 380);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCAU"), &NucleotideVec::from_lossy("AC")), 330);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CAAG")), 280);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AU"), &NucleotideVec::from_lossy("AAAU")), 380);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GAAU")), 330);

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CCAG")), 280);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AU"), &NucleotideVec::from_lossy("ACAU")), 380);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GCAU")), 330);
    }

    #[test]
    fn test_vrna_interior_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACA"), &NucleotideVec::from_lossy("UGAAU")), 370);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACAA"), &NucleotideVec::from_lossy("UGAAU")), 290);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GUAGU"), &NucleotideVec::from_lossy("AGGC")), 260);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AUAGU"), &NucleotideVec::from_lossy("AGGU")), 330);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GGC"), &NucleotideVec::from_lossy("GUGC")), 110);
    }


    #[test]
    fn test_vrna_bulge_n_evaluation() {
        let model = ViennaRNA::default();

        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAAAAAG"), &NucleotideVec::from_lossy("CG")), 440);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAAAAAAAG"), &NucleotideVec::from_lossy("CG")), 470);
    }

    #[test]
    fn test_vrna_multibranch() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AAAA");
        let seg2 = &NucleotideVec::from_lossy("AAA");
        let seg3 = &NucleotideVec::from_lossy("AAAAA");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 720);

        let seg1 = &NucleotideVec::from_lossy("GAAC");
        let seg2 = &NucleotideVec::from_lossy("GAC");
        let seg3 = &NucleotideVec::from_lossy("GAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 330);

        let seg1 = &NucleotideVec::from_lossy("GAAC");
        let seg2 = &NucleotideVec::from_lossy("GAC");
        let seg3 = &NucleotideVec::from_lossy("GAAAAAAAAAAAAAAAAAAC");
        let energy = model.multibranch(&[seg1, seg2, seg3]);
        assert_eq!(energy, 330);

    }

    #[test]
    fn test_vrna_exterior_single_branch() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::from_lossy("UG");
        let seg2 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120); 

        let seg1 = &NucleotideVec::from_lossy("G");
        let seg2 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, -120);
 
        let seg1 = &NucleotideVec::from_lossy("UG");
        let seg2 = &NucleotideVec::from_lossy("C");
        let energy = model.exterior(&[seg1, seg2]);
        assert_eq!(energy, 0); 
    }

    #[test]
    fn test_vrna_exterior_two_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -240);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUUG");
        let seg3 = &NucleotideVec::from_lossy("C");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -120);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CG");
        let seg3 = &NucleotideVec::from_lossy("CU");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -290);

        let seg1 = &NucleotideVec::from_lossy("ACA");
        let seg2 = &NucleotideVec::from_lossy("UGG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3]);
        assert_eq!(energy, -130);
    }

    #[test]
    fn test_vrna_exterior_three_branches() {
        let model = ViennaRNA::default();

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("CUG");
        let seg4 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]);
        assert_eq!(energy, -360);

        let seg1 = &NucleotideVec::from_lossy("AUG");
        let seg2 = &NucleotideVec::from_lossy("CUG");
        let seg3 = &NucleotideVec::from_lossy("UUG");
        let seg4 = &NucleotideVec::from_lossy("CUG");
        let energy = model.exterior(&[seg1, seg2, seg3, seg4]);
        assert_eq!(energy, -240);
    }

 
    #[test]
    fn test_evaluations() {
        let model = ViennaRNA::default();

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "GAAAAG";
        let dbr = "(....)";
        let e37 = 630;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "ACGUUAAAGACGU";
        let dbr = "(((((...)))))";
        let e37 = -170;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "AGACGACAAGGUUGAAUCGC";
        let dbr = ".(.(((.(....)...))))";
        let e37 = 420;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "GAGUAGUGGAACCAGGCUAU";
        let dbr = ".((...((....))..))..";
        let e37 = 190;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);

        let seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACC";
        let dbr = "...........(.(((((........)))))..)......";
        let e37 = -210;
        assert_eq!(model.energy_of_structure(&NucleotideVec::from_lossy(seq), &PairTable::try_from(dbr).expect("valid")), e37);



    }
}

 
