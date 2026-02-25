
use std::fmt;
use std::borrow::Borrow;
use std::ops::Deref;

use log::warn;
use colored::*;


#[derive(Debug)]
pub enum SequenceError {
    Plain(String),
    InvalidChar(char),
    Separator(char),
}

impl fmt::Display for SequenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SequenceError::Plain(s) => {
                write!(f, "ERROR: {}", s)
            }
            SequenceError::InvalidChar(c) => {
                write!(f, "Unsupported nucleotide: '{}'", c)
            }
            SequenceError::Separator(c) => {
                write!(f, "Unexpected strand separation character '{}'", c)
            }
        }
    }
}

impl std::error::Error for SequenceError {}


#[repr(u8)]
#[derive(Clone, Hash, Copy, Debug, Eq, PartialEq)]
pub enum Base { A, C, G, U, N, StrandBreak }
pub const BCOUNT: usize = 5; // 5 Base variants for tables.

impl TryFrom<char> for Base {
    type Error = SequenceError;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'U' | 'T' => Ok(Base::U),
            'N' => Ok(Base::N),
            '&' | '+' => Ok(Base::StrandBreak),
            _ => Err(SequenceError::InvalidChar(c)),
        }
    }
}

impl fmt::Display for Base {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let c = match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::U => 'U',
            Base::N => 'N',
            Base::StrandBreak => '+',
        };
        write!(f, "{}", c)
    }
}


#[derive(Clone, Hash, Debug, Eq, PartialEq)]
pub struct NucleotideVec(pub Vec<Base>);

impl Deref for NucleotideVec {
    type Target = [Base];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl Borrow<[Base]> for NucleotideVec {
    fn borrow(&self) -> &[Base] {
        &self.0
    }
}

impl TryFrom<&str> for NucleotideVec {
    type Error = SequenceError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut vec = Vec::with_capacity(s.len());
        for c in s.chars() {
            vec.push(Base::try_from(c)?);
        }
        Ok(NucleotideVec(vec))
    }
}

impl fmt::Display for NucleotideVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for base in &self.0 {
            write!(f, "{}", base)?;
        }
        Ok(())
    }
}

impl NucleotideVec {
    pub fn from_lossy(s: &str) -> Self {
        let vec = s.chars().map(|c| {
            Base::try_from(c).unwrap_or_else(|e| {
                warn!("{} {} -> converted to 'N'", "WARNING:".red(), e);
                Base::N
            })
        }).collect();
        NucleotideVec(vec)
    }

    pub fn has_indistinguishable_strands(&self) -> bool {
        let blocks: Vec<&[Base]> = self.0
            .split(|b| *b == Base::StrandBreak)
            .collect();
        blocks.iter().enumerate().any(|(i, a)| {
            blocks.iter().skip(i + 1).any(|b| a == b)
        })
    }

}

const PAIR_LOOKUP: [[PairTypeRNA; BCOUNT]; BCOUNT] = {
    use Base::*;
    use PairTypeRNA::*;
    let mut table = [[NN; BCOUNT]; BCOUNT];
    table[A as usize][U as usize] = AU;
    table[U as usize][A as usize] = UA;
    table[C as usize][G as usize] = CG;
    table[G as usize][C as usize] = GC;
    table[G as usize][U as usize] = GU;
    table[U as usize][G as usize] = UG;
    table
};

#[repr(usize)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum PairTypeRNA {
    AU = 0,
    UA = 1,
    CG = 2,
    GC = 3,
    GU = 4,
    UG = 5,
    NN = 6,
}
pub const PCOUNT: usize = 7; // 7 Pair variants for tables.

impl From<(Base, Base)> for PairTypeRNA {
    fn from(pair: (Base, Base)) -> Self {
        PAIR_LOOKUP[pair.0 as usize][pair.1 as usize]
    }
}

impl fmt::Display for PairTypeRNA {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            PairTypeRNA::AU => "A-U",
            PairTypeRNA::UA => "U-A",
            PairTypeRNA::CG => "C-G",
            PairTypeRNA::GC => "G-C",
            PairTypeRNA::GU => "G-U",
            PairTypeRNA::UG => "U-G",
            PairTypeRNA::NN => "N-N",
        };
        write!(f, "{}", s)
    }
}

impl PairTypeRNA {
    pub fn new(pair: (Base, Base)) -> Self {
        let pt = PAIR_LOOKUP[pair.0 as usize][pair.1 as usize];

        if pt == PairTypeRNA::NN {
            warn!("{} Invalid base pair: {}-{} -> converted to {}", "WARNING:".red(), pair.0, pair.1, pt);
        }

        pt
    }
    
    pub fn is_ru(&self) -> bool {
       matches!(self
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA)
    }

    pub fn is_wcf(&self) -> bool {
       matches!(self
            , PairTypeRNA::GC | PairTypeRNA::CG 
            | PairTypeRNA::AU | PairTypeRNA::UA)
    }

    pub fn is_wobble(&self) -> bool {
       matches!(self, PairTypeRNA::GU | PairTypeRNA::UG)
    }

    pub fn can_pair(&self) -> bool {
       self != &PairTypeRNA::NN
    }
    
    pub fn invert(&self) -> PairTypeRNA {
        use PairTypeRNA::*;
        match self {
            AU => UA,
            UA => AU,
            CG => GC,
            GC => CG,
            GU => UG,
            UG => GU,
            NN => NN,
        }
    }
}


