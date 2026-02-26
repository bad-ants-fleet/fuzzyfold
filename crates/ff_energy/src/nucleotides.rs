use std::fmt;
use std::borrow::Borrow;
use std::ops::Deref;

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


/// All currently supported Bases. 
///
/// Note the following choices:
///     - we distinguish U, T, and, PU.
///     - strand break (SB) is treated as a Base, but it is not part of BCOUNT,
///     because input sequences are parsed into slices of Bases (NucleotideVec)
///     however, they are always removed during energy evaluation, before table
///     lookups, so it's not part of BCOUNT.
///     - 
#[repr(u8)]
#[derive(Clone, Hash, Copy, Debug, Eq, PartialEq)]
pub enum Base { A, C, G, U, T, PU, SB }
pub const BCOUNT: usize = 6; // all except SB

impl Base {
    /// T is always assumed to be U (for now).
    /// PU always uses U parameters if unpaired.
    #[inline(always)]
    pub const fn canonical_rna_index(self) -> usize {
        match self {
            Base::T => Base::U as usize,
            Base::PU => Base::U as usize,
            _ => self as usize,
        }
    }
}

impl TryFrom<char> for Base {
    type Error = SequenceError;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c.to_ascii_uppercase() {
            'A' => Ok(Base::A),
            'C' => Ok(Base::C),
            'G' => Ok(Base::G),
            'U' => Ok(Base::U),
            'T' => Ok(Base::T),
            'P' => Ok(Base::PU),
            '&' | '+' => Ok(Base::SB),
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
            Base::T => 'T',
            Base::PU => 'P',
            Base::SB => '+',
        };
        write!(f, "{}", c)
    }
}

/// The main representation of a (multi-stranded) nucleic acid sequence.
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
            Base::try_from(c).unwrap_or_else(|c| panic!("There's an unknown character {} in your sequence.", c))
        }).collect();
        NucleotideVec(vec)
    }

    pub fn has_indistinguishable_strands(&self) -> bool {
        let blocks: Vec<&[Base]> = self.0
            .split(|b| *b == Base::SB)
            .collect();
        blocks.iter().enumerate().any(|(i, a)| {
            blocks.iter().skip(i + 1).any(|b| a == b)
        })
    }

}

/// Now we are in RNA territory. Stacking tables now distinguish AP PA pairs
/// from AU UA pairs, that's why they are listed explcitly 
/// (in contrast to AG / GA pairs which are treated always as UG / GU.
#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum PairTypeRNA {
    AU = 0,
    UA = 1,
    CG = 2,
    GC = 3,
    GU = 4,
    UG = 5,
    AP = 6,
    PA = 7,
    NN = 8,
}
pub const PCOUNT: usize = 8; // all except SB

const PAIR_LOOKUP: [[PairTypeRNA; BCOUNT]; BCOUNT] = {
    use Base::*;
    use PairTypeRNA::*;
    let mut table = [[NN; BCOUNT]; BCOUNT];
    table[A as usize][U as usize] = AU;
    table[U as usize][A as usize] = UA;
    table[A as usize][T as usize] = AU;
    table[T as usize][A as usize] = UA;
    table[C as usize][G as usize] = CG;
    table[G as usize][C as usize] = GC;
    table[G as usize][U as usize] = GU;
    table[U as usize][G as usize] = UG;
    table[G as usize][T as usize] = GU;
    table[T as usize][G as usize] = UG;
    table[A as usize][PU as usize] = AP;
    table[PU as usize][A as usize] = PA;
    table[G as usize][PU as usize] = GU;
    table[PU as usize][G as usize] = UG;
    table
};

const FALLBACK_LOOKUP: [[PairTypeRNA; BCOUNT]; BCOUNT] = {
    use Base::*;
    use PairTypeRNA::*;
    let mut table = [[NN; BCOUNT]; BCOUNT];
    table[A as usize][U as usize] = AU;
    table[U as usize][A as usize] = UA;
    table[A as usize][T as usize] = AU;
    table[T as usize][A as usize] = UA;
    table[C as usize][G as usize] = CG;
    table[G as usize][C as usize] = GC;
    table[G as usize][U as usize] = GU;
    table[U as usize][G as usize] = UG;
    table[G as usize][T as usize] = GU;
    table[T as usize][G as usize] = UG;
    table[A as usize][PU as usize] = AU;
    table[PU as usize][A as usize] = UA;
    table[G as usize][PU as usize] = GU;
    table[PU as usize][G as usize] = UG;
    table
};

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
            PairTypeRNA::AP => "A-P",
            PairTypeRNA::PA => "P-A",
            PairTypeRNA::NN => "N-N",
        };
        write!(f, "{}", s)
    }
}

const PAIR_INVERT: [PairTypeRNA; PCOUNT + 1] = {
    use PairTypeRNA::*;
    [UA, AU, GC, CG, UG, GU, PA, AP, NN]
};

impl PairTypeRNA {
    pub fn from_fallback(pair: (Base, Base)) -> Self {
        FALLBACK_LOOKUP[pair.0 as usize][pair.1 as usize]
    }

    pub fn is_ru(&self) -> bool {
       matches!(self
            , PairTypeRNA::GU | PairTypeRNA::UG 
            | PairTypeRNA::AU | PairTypeRNA::UA)
    }

    /// TODO: Is not used at the moment.
    pub fn is_ap(&self) -> bool {
       matches!(self, PairTypeRNA::AP | PairTypeRNA::PA)
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
        PAIR_INVERT[*self as usize]
    }
}


