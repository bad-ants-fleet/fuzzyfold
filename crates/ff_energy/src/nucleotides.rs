//! Nucleotides and base pairs: `Base`, `NucleotideVec` and `PairTypeRNA`.

use std::fmt;
use std::borrow::Borrow;
use std::ops::Deref;

/// Errors when parsing a sequence.
#[derive(Debug)]
pub enum SequenceError {
    /// A character that is not a nucleotide or a strand break.
    InvalidBase(char),
    /// A nucleotide that is not allowed in RNA (T).
    InvalidRnaBase(Base),
    /// A nucleotide that is not allowed in DNA (U, P, I).
    InvalidDnaBase(Base),
    /// A strand without bases, ended at this position.
    EmptyStrand(usize),
}

impl fmt::Display for SequenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SequenceError::InvalidBase(c) => {
                write!(f, "Unsupported nucleotide: '{}'", c)
            }
            SequenceError::InvalidRnaBase(c) => {
                write!(f, "Unsupported nucleotide for RNA: '{}'", c)
            }
            SequenceError::InvalidDnaBase(c) => {
                write!(f, "Unsupported nucleotide for DNA: '{}'", c)
            }
            SequenceError::EmptyStrand(i) => {
                write!(f, "Empty strand at position {}", i)
            }
        }
    }
}

impl std::error::Error for SequenceError {}


/// A nucleotide, or the end of a strand.
///
/// Modified bases have two-letter names (`PU`, `IN`). U and T are distinct
/// bases, but share their parameters. A strand break
/// or strand wrap takes a position in a sequence, like `+` or `^` in a
/// structure, so that indices of sequence and structure agree. It never
/// pairs and is not looked up in parameter tables.
#[repr(u8)]
#[derive(Clone, Hash, Copy, Debug, Eq, PartialEq)]
pub enum Base {
    A,
    C,
    G,
    U,
    T,
    /// Pseudouridine, written `P`.
    PU,
    /// Inosine, written `I`.
    IN,
    /// A strand break, written `+` (also parsed from `&`): the strand before
    /// it ends.
    SB,
    /// A strand wrap, written `^`: the strand before it is circular, it
    /// continues at its own first base.
    SW,
}

/// The number of bases in parameter tables: all except `SB` and `SW`.
const BCOUNT: usize = 7;

impl Base {
    /// The index of the base in tables for unpaired bases (dangles,
    /// mismatches): T and P use the values of U, I those of G.
    #[inline(always)]
    pub(crate) const fn canonical_rna_index(self) -> usize {
        match self {
            Base::T => Base::U as usize,
            Base::PU => Base::U as usize,
            Base::IN => Base::G as usize,
            _ => self as usize,
        }
    }
}

/// Case insensitive.
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
            'I' => Ok(Base::IN),
            '&' | '+' => Ok(Base::SB),
            '^' => Ok(Base::SW),
            _ => Err(SequenceError::InvalidBase(c)),
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
            Base::IN => 'I',
            Base::SB => '+',
            Base::SW => '^',
        };
        write!(f, "{}", c)
    }
}

/// A sequence of one or more strands. A strand ends at a strand break
/// (`Base::SB`), at a strand wrap (`Base::SW`, a circular strand), or at
/// the end of the sequence.
///
/// `TryFrom<&str>` accepts all bases; `try_from_rna` and `try_from_dna` also
/// check that the bases fit the nucleic acid. Every strand has at least one
/// base: a strand break or strand wrap must follow one, and a sequence
/// must not end with a strand break.
///
/// ```
/// use ff_energy::{Base, NucleotideVec};
///
/// let seq = NucleotideVec::try_from_rna("GGGA+UCCC").unwrap();
/// assert_eq!(seq[4], Base::SB);
/// assert!(NucleotideVec::try_from_rna("GGTA").is_err());
/// ```
#[derive(Clone, Hash, Debug, Eq, PartialEq)]
pub struct NucleotideVec(Vec<Base>);

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
        NucleotideVec::checked(vec)
    }
}

/// A sequence from bases, e.g. a modified copy of another sequence; every
/// strand must have at least one base.
impl TryFrom<Vec<Base>> for NucleotideVec {
    type Error = SequenceError;

    fn try_from(vec: Vec<Base>) -> Result<Self, Self::Error> {
        NucleotideVec::checked(vec)
    }
}

/// The number of circular strands of a sequence.
pub(crate) fn circular_strands(sequence: &[Base]) -> usize {
    sequence.iter().filter(|&&b| b == Base::SW).count()
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
    /// An RNA sequence: T is an error.
    pub fn try_from_rna(s: &str) -> Result<Self, SequenceError> {
        let vec = s.chars()
            .map(|c| {
                let base = Base::try_from(c)?;
                match base {
                    Base::T => Err(SequenceError::InvalidRnaBase(base)),
                    _ => Ok(base),
                }
            })
        .collect::<Result<Vec<_>, _>>()?;
        NucleotideVec::checked(vec)
    }

    /// A DNA sequence: U, P and I are errors.
    pub fn try_from_dna(s: &str) -> Result<Self, SequenceError> {
        let vec = s.chars()
            .map(|c| {
                let base = Base::try_from(c)?;
                match base {
                    Base::U => Err(SequenceError::InvalidDnaBase(base)),
                    Base::PU | Base::IN => Err(SequenceError::InvalidDnaBase(base)),
                    _ => Ok(base),
                }
            })
        .collect::<Result<Vec<_>, _>>()?;
        NucleotideVec::checked(vec)
    }

    /// The sequence, if every strand has at least one base.
    fn checked(vec: Vec<Base>) -> Result<Self, SequenceError> {
        for (i, &b) in vec.iter().enumerate() {
            if !matches!(b, Base::SB | Base::SW) {
                continue;
            }
            let empty = i == 0 || matches!(vec[i - 1], Base::SB | Base::SW);
            if empty || (b == Base::SB && i + 1 == vec.len()) {
                return Err(SequenceError::EmptyStrand(i));
            }
        }
        Ok(NucleotideVec(vec))
    }

    /// True if at least one strand is circular.
    pub fn has_circular_strand(&self) -> bool {
        self.0.contains(&Base::SW)
    }

    /// True if two strands have the same sequence, and both or none are
    /// circular. Their order cannot be told apart, which matters for
    /// symmetry corrections.
    pub fn has_indistinguishable_strands(&self) -> bool {
        // Each strand with its strand wrap, if it is circular.
        let blocks: Vec<&[Base]> = self.0
            .split_inclusive(|b| matches!(b, Base::SB | Base::SW))
            .map(|s| s.strip_suffix(&[Base::SB]).unwrap_or(s))
            .collect();
        blocks.iter().enumerate().any(|(i, a)| {
            blocks.iter().skip(i + 1).any(|b| a == b)
        })
    }

}

/// The type of a base pair, 5' base first, or `NN` if the bases cannot pair.
///
/// Pairs with a modified base have their own types where the parameters
/// distinguish them (A-P, and the pairs of inosine). G-P is a G-U pair, and
/// T pairs like U.
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
    IU = 8,
    UI = 9,
    CI = 10,
    IC = 11,
    NN = 12,
}

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
    table[IN as usize][U as usize] = IU;
    table[U as usize][IN as usize] = UI;
    table[C as usize][IN as usize] = CI;
    table[IN as usize][C as usize] = IC;
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
    table[IN as usize][U as usize] = GU;
    table[U as usize][IN as usize] = UG;
    table[C as usize][IN as usize] = CG;
    table[IN as usize][C as usize] = GC;
    table
};

/// Panics if a base is `Base::SB` or `Base::SW`.
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
            PairTypeRNA::IU => "I-U",
            PairTypeRNA::UI => "U-I",
            PairTypeRNA::CI => "C-I",
            PairTypeRNA::IC => "I-C",
            PairTypeRNA::NN => "N-N",
        };
        write!(f, "{}", s)
    }
}

impl PairTypeRNA {
    /// The pair type with modified bases replaced by the canonical ones
    /// (A-P as A-U, inosine as G), for parameters that do not distinguish
    /// them: hairpins, interior loops, mismatches.
    pub(crate) fn from_fallback(pair: (Base, Base)) -> Self {
        FALLBACK_LOOKUP[pair.0 as usize][pair.1 as usize]
    }

    /// False for `NN`.
    pub fn can_pair(&self) -> bool {
       self != &PairTypeRNA::NN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_wrap() {
        let seq = NucleotideVec::try_from("GGGA^UCCC").unwrap();
        assert_eq!(seq[4], Base::SW);
        assert_eq!(seq.to_string(), "GGGA^UCCC");
        assert!(seq.has_circular_strand());
        for s in ["^GGG", "GGG+^CC", "GG^^", "+GGG", "GG++CC", "GGG+"] {
            assert!(matches!(NucleotideVec::try_from(s), Err(SequenceError::EmptyStrand(_))), "{s}");
        }
        assert!(NucleotideVec::try_from_rna("^GGG").is_err());
        assert!(NucleotideVec::try_from_dna("^GGG").is_err());
    }

    #[test]
    fn test_indistinguishable_strands() {
        let strands = |s: &str| NucleotideVec::try_from(s).unwrap().has_indistinguishable_strands();
        assert!(strands("GGC+GGC"));
        assert!(strands("GGC^GGC^"));
        assert!(!strands("GGC^GGC"));
        assert!(!strands("GGC+GCC"));
    }

    #[test]
    fn test_dna_bases() {
        assert!(NucleotideVec::try_from_dna("ACGT").is_ok());
        for s in ["ACGU", "ACGP", "ACGI"] {
            assert!(matches!(NucleotideVec::try_from_dna(s), Err(SequenceError::InvalidDnaBase(_))), "{s}");
        }
    }
}
