use paste::paste;
use log::info;
use colored::*; 
use crate::EnergyModel;
use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;
use crate::parameters::*;
use crate::PairTypeRNA;
use crate::K0;
const T_REF: f64 = 37.0;

pub enum TableSource<T: 'static> {
    Static(&'static T),
    Owned(Box<T>),
}

impl<T> TableSource<T> {
    #[inline]
    pub fn get(&self) -> &T {
        match self {
            TableSource::Static(t) => t,
            TableSource::Owned(b) => b,
        }
    }
}

pub struct ViennaRNA {
    min_hp_size: usize,
    temperature: f64,
    params: &'static ThermoParams,

    stack: StackParams,
    mismatch_hairpin: MismatchParams,
    mismatch_interior: MismatchParams,
    mismatch_interior_1n: MismatchParams,
    mismatch_interior_23: MismatchParams,
    mismatch_multi: MismatchParams,
    mismatch_exterior: MismatchParams,
                                   
    dangle5: DangleParams,
    dangle3: DangleParams,
                                                                              
    int11: TableSource<Int11Params>,
    int21: TableSource<Int21Params>,
    int22: TableSource<Int22Params>,

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

macro_rules! rescale_table {
    ($self:ident, $field:ident, $scale:ident) => {{
        paste! {
            $self.params.[<$field _en37>]
                .rescale($self.params.[<$field _enth>], $scale)
        }
    }};
}

macro_rules! rescale_value {
    ($self:ident, $field:ident, $scale:ident) => {
        paste! {
            $self.$field =
                $self.params.[<$field _en37>]
                    .rescale(&$self.params.[<$field _enth>], $scale);
        }
    };
}

macro_rules! copy_table {
    ($self:ident, $field:ident) => {
        paste! {
            $self.$field = *$self.params.[<$field _en37>];
        }
    };
}

macro_rules! copy_value {
    ($self:ident, $field:ident) => {
        paste! {
            $self.$field = $self.params.[<$field _en37>];
        }
    };
}

impl Default for ViennaRNA {
    fn default() -> Self {
        Self::new(&TURNER2004)
    }
}

impl ViennaRNA {
    pub fn new(params: &'static ThermoParams) -> Self {
        Self {
            min_hp_size: 3,
            temperature: 37.0,
            params,

            stack: *params.stack_en37,
            mismatch_hairpin: *params.mismatch_hairpin_en37,
            mismatch_interior: *params.mismatch_interior_en37,
            mismatch_interior_1n: *params.mismatch_interior_1n_en37,
            mismatch_interior_23: *params.mismatch_interior_23_en37,
            mismatch_multi: *params.mismatch_multi_en37,
            mismatch_exterior: *params.mismatch_exterior_en37,

            dangle5: *params.dangle5_en37,
            dangle3: *params.dangle3_en37,

            int11: TableSource::Static(params.int11_en37),
            int21: TableSource::Static(params.int21_en37),
            int22: TableSource::Static(params.int22_en37),

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

            triloops: params.triloops.to_vec(),
            tetraloops: params.tetraloops.to_vec(),
            hexaloops: params.hexaloops.to_vec(),
        }
    }

    // NOTE: Setting temperature resets the tables!
    pub fn reset_with_temperature(&mut self, temp: f64) {
        self.temperature = temp;
        if (temp - T_REF).abs() < 1e-4 {
            copy_table!(self, stack);
            copy_table!(self, mismatch_hairpin);
            copy_table!(self, mismatch_interior);
            copy_table!(self, mismatch_interior_1n);
            copy_table!(self, mismatch_interior_23);
            copy_table!(self, mismatch_multi);
            copy_table!(self, mismatch_exterior);
            copy_table!(self, dangle5);
            copy_table!(self, dangle3);
            self.int11 = TableSource::Static(self.params.int11_en37);
            self.int21 = TableSource::Static(self.params.int21_en37);
            self.int22 = TableSource::Static(self.params.int22_en37);
            copy_table!(self, hairpin);
            copy_table!(self, bulge);
            copy_table!(self, interior);

            copy_value!(self, duplex_init);
            copy_value!(self, terminal_ru);
            self.lxc = self.params.lxc;
            copy_value!(self, ninio);
            self.ninio_max = self.params.ninio_max;
            copy_value!(self, ml_base);
            copy_value!(self, ml_closing);
            copy_value!(self, ml_intern);

            self.triloops = self.params.triloops.to_vec();
            self.tetraloops = self.params.tetraloops.to_vec();
            self.hexaloops = self.params.hexaloops.to_vec();
            return;
        }

        let kelvin = temp + K0;
        let scale = kelvin / (T_REF + K0);
        rescale_table!(self, stack, scale);
        rescale_table!(self, mismatch_hairpin, scale);
        rescale_table!(self, mismatch_interior, scale);
        rescale_table!(self, mismatch_interior_1n, scale);
        rescale_table!(self, mismatch_interior_23, scale);
        rescale_table!(self, mismatch_multi, scale);
        rescale_table!(self, mismatch_exterior, scale);
        rescale_table!(self, dangle5, scale);
        rescale_table!(self, dangle3, scale);
        self.int11 = TableSource::Owned(Box::new(rescale_table!(self, int11, scale)));
        self.int21 = TableSource::Owned(Box::new(rescale_table!(self, int21, scale)));
        self.int22 = TableSource::Owned(Box::new(rescale_table!(self, int22, scale)));
        rescale_table!(self, hairpin, scale);
        rescale_table!(self, bulge, scale);
        rescale_table!(self, interior, scale);

        rescale_value!(self, duplex_init, scale);
        rescale_value!(self, terminal_ru, scale);
        self.lxc = self.params.lxc * scale;
        rescale_value!(self, ninio, scale);
        self.ninio_max = self.params.ninio_max; // NOTE
        rescale_value!(self, ml_base, scale);
        rescale_value!(self, ml_closing, scale);
        rescale_value!(self, ml_intern, scale);

        self.triloops = self.params.triloops
            .iter().map(|e| e.rescaled(scale)).collect();
        self.tetraloops = self.params.tetraloops
            .iter().map(|e| e.rescaled(scale)).collect();
        self.hexaloops = self.params.hexaloops
            .iter().map(|e| e.rescaled(scale)).collect();
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
            .find(|e| {
                let bytes = e.seq.as_bytes();
                for (b, &c) in seq.iter().zip(bytes.iter()) {
                    let base = match c {
                        b'A' => Base::A,
                        b'C' => Base::C,
                        b'G' => Base::G,
                        b'U' => Base::U,
                        _ => return false,
                    };

                    if *b != base {
                        return false;
                    }
                }
                true
            })
            .map(|e| e.g37)
    }

    // TODO: Return result!
    fn hairpin(&self, seq: &[Base]) -> i32 {
        let n = seq.len() - 2;
        if n < self.min_hp_size {
            panic!("Invalid hairpin size {n}");
        }

        // Special hairpin energies
        if let Some(en) = self.hairpin_bonus(seq) {
            return en;
        }

        let closing = PairTypeRNA::new((seq[0], *seq.last().unwrap()));
        assert!(closing.can_pair());

        // Initiation terms
        let mut en = if n <= 30 {
            self.hairpin[n]
        } else {
            self.hairpin[30] + (self.lxc * ((n as f64) / 30.).ln()) as i32
        };

        if n == 3 && closing.is_ru() {
            en += self.terminal_ru;
        } else if n > 3 {
            //println!("mmh {} {} {}", closing, seq[1], seq[n]);
            en += self.mismatch_hairpin
                [closing as usize]
                [seq[1] as usize]
                [seq[n] as usize];
        }
        en
    }

    fn interior(&self, fwdseq: &[Base], revseq: &[Base]) -> i32 {
        let outer = PairTypeRNA::new((*fwdseq.first().unwrap(), *revseq.last().unwrap()));
        let inner = PairTypeRNA::from((*revseq.first().unwrap(), *fwdseq.last().unwrap()));
        assert!(outer.can_pair() && inner.can_pair());

        match (fwdseq.len(), revseq.len()) {
            (2, 2) => 
                self.stack[outer as usize][inner as usize],
            (3, 2) | (2, 3) => { //NOTE: SpecialC if C adjacent to paired C missing!
                self.bulge[1] + 
                self.stack[outer as usize][inner as usize]},
            (3, 3) => 
                self.int11.get()[outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize],
            (3, 4) => 
                self.int21.get()
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][revseq[1] as usize]
                [revseq[2] as usize],
            (4, 3) => 
                self.int21.get()
                [inner as usize][outer as usize]
                [revseq[1] as usize][fwdseq[1] as usize]
                [fwdseq[2] as usize],
            (4, 4) => 
                self.int22.get()
                [outer as usize][inner as usize]
                [fwdseq[1] as usize][fwdseq[2] as usize]
                [revseq[1] as usize][revseq[2] as usize],
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
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[revseq.len() - 2] as usize] +
                    self.mismatch_interior_1n
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[fwdseq.len() - 2] as usize];

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
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[revseq.len() - 2] as usize] +
                    self.mismatch_interior_23
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[fwdseq.len() - 2] as usize];
                en += self.ninio;
                en += self.interior[5];
                en
            }
            (lfwd, lrev) => { 
                let mut en = self.mismatch_interior
                    [outer as usize][fwdseq[1] as usize]
                    [revseq[lrev - 2] as usize] +
                    self.mismatch_interior
                    [inner as usize][revseq[1] as usize]
                    [fwdseq[lfwd - 2] as usize];

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
        }
    }

    fn multibranch(&self, segments: &[&[Base]]) -> i32 {
        // For warning purposes only.
        let closing = PairTypeRNA::new((segments[0][0], *segments.last().unwrap().last().unwrap()));
        assert!(closing.can_pair());

        // Number of stems in the multiloop.
        let n = segments.len(); 

        let mut en = 0;
        for i in 0..n {
            let j = (i + 1) % n; 
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            assert!(pair.can_pair());
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
                    [pair as usize][b5 as usize][b3 as usize],
                (Some(&b5), None) => 
                    self.dangle5
                     [pair as usize][b5 as usize].min(0),
                (None, Some(&b3)) => 
                    self.dangle3
                    [pair as usize][b3 as usize].min(0),
                _ => 0,
            };
            en += den;
        }
 
        // Number of unpaired bases in the multiloop.
        let m: usize = segments.iter().map(|s| s.len() - 2).sum();
        en + self.ml_base * m as i32
           + self.ml_closing
           + self.ml_intern * n as i32
    }

    fn exterior(&self, segments: &[&[Base]]) -> i32 {
        let mut en = 0;
        let n = segments.len() - 1; 
        for i in 0..n {
            let j = i + 1; 
            
            let pair = PairTypeRNA::from((*segments[i].last().unwrap(), segments[j][0]));
            assert!(pair.can_pair());
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
                    [pair as usize][b5 as usize][b3 as usize],
                (Some(&b5), None) => 
                    self.dangle5
                    [pair as usize][b5 as usize].min(0),
                (None, Some(&b3)) => 
                     self.dangle3
                    [pair as usize][b3 as usize].min(0),
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
                self.exterior(&slices) + self.duplex_init
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
    use crate::parameters::TURNER2004;
    use crate::NucleotideVec;

    #[test]
    fn test_vrna_hairpin_evaluation() {
        let model = ViennaRNA::new(&TURNER2004);
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
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAU")), 530);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAAU")), 580);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("ACCCCU")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("ACCCCCU")), 550);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAU")), 540);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAU")), 510);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAAAAU")), 610);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(30)))), 620);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(31)))), 623);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(32)))), 626);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(33)))), 630);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(34)))), 633);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy(&format!("C{}G", "A".repeat(35)))), 636);
    }

    fn test_cannot_pair() {
        let model = ViennaRNA::new(&TURNER2004);
        //assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAG")), 590);
        //assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAAC")), 590);
        //assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAAAAA")), 660);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("GAAAG")), 9999);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("CAAAC")), 9999);
        assert_eq!(model.hairpin(&NucleotideVec::from_lossy("AAAAAAAAAAA")), 9999);
    }

    #[test]
    fn test_vrna_stacking_evaluation() {
        let model = ViennaRNA::new(&TURNER2004);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CG"), &NucleotideVec::from_lossy("CG")), -240);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AC"), &NucleotideVec::from_lossy("GU")), -220);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GU"), &NucleotideVec::from_lossy("AC")), -220);
    }

    #[test]
    fn test_vrna_int11_evaluation() {
        let model = ViennaRNA::new(&TURNER2004);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CCG"), &NucleotideVec::from_lossy("CGG")), 50);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAG"), &NucleotideVec::from_lossy("CAG")), 90);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACU"), &NucleotideVec::from_lossy("AAU")), 190);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCU"), &NucleotideVec::from_lossy("AUC")), 120);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GCU"), &NucleotideVec::from_lossy("AGC")), 120);
    }

    #[test]
    fn test_vrna_int21_evaluation() {
        let model = ViennaRNA::new(&TURNER2004);
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
        let model = ViennaRNA::new(&TURNER2004);
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
        let model = ViennaRNA::new(&TURNER2004);
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
        let model = ViennaRNA::new(&TURNER2004);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACA"), &NucleotideVec::from_lossy("UGAAU")), 370);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("ACAA"), &NucleotideVec::from_lossy("UGAAU")), 290);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GUAGU"), &NucleotideVec::from_lossy("AGGC")), 260);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("AUAGU"), &NucleotideVec::from_lossy("AGGU")), 330);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("GGC"), &NucleotideVec::from_lossy("GUGC")), 110);
    }

    #[test]
    fn test_vrna_bulge_n_evaluation() {
        let model = ViennaRNA::new(&TURNER2004);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAAAAAG"), &NucleotideVec::from_lossy("CG")), 440);
        assert_eq!(model.interior(&NucleotideVec::from_lossy("CAAAAAAAAG"), &NucleotideVec::from_lossy("CG")), 470);
    }

    #[test]
    fn test_vrna_multibranch() {
        let model = ViennaRNA::new(&TURNER2004);
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
        let model = ViennaRNA::new(&TURNER2004);

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
        let model = ViennaRNA::new(&TURNER2004);

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
        let model = ViennaRNA::new(&TURNER2004);

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
        let model = ViennaRNA::new(&TURNER2004);

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
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
 
    #[test]
    fn test_multi_evaluations() {
        let model = ViennaRNA::new(&TURNER2004);

        let seq = "GAAAAC";
        let dbr = "(....)";
        let e37 = 450;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")), e37);

        let fseq = "GA+AAAC";
        let fdbr = "(.+...)";
 
        let rseq = "AAAC+GA";
        let rdbr = "...(+).";
        assert_eq!(
            model.energy_of_structure(
                &NucleotideVec::from_lossy(fseq), 
                &MultiPairTable::try_from(fdbr).expect("valid")), 
            model.energy_of_structure(
                &NucleotideVec::from_lossy(rseq), 
                &MultiPairTable::try_from(rdbr).expect("valid"))
        );

        let seq = "GAA+AAC";
        let dbr = "(..+..)";
        let e37 = 300;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")), e37);

        let seq = "GC+UUUUAGU+AU+AC";
        let dbr = "((+(...)).+..+.)";
        let e37 = 1140;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")), e37);

        let seq = "GC&UUUUAGU&AGAAACU&AGAAACU&AC";
        let dbr = "((&(...)).&.(...).&.(...).&.)";
        let e37 = 2020;
        assert_eq!(model.energy_of_structure(
                &NucleotideVec::from_lossy(seq), 
                &MultiPairTable::try_from(dbr).expect("valid")), e37);
 
    }

}

