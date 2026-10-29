//! Dot-bracket notation with pseudoknot support.

use std::fmt;
use std::ops::Deref;
use std::ops::DerefMut;
use std::convert::TryFrom;
use crate::error::StructureError;
use crate::MultiPairTable;
use crate::MultiStruct;

/// The opening and closing bracket of each level.
const OPENING: &[u8; DotBracket::LEVELS] = b"([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const CLOSING: &[u8; DotBracket::LEVELS] = b")]}>abcdefghijklmnopqrstuvwxyz";

/// One character of dot-bracket notation, stored as its ASCII byte.
///
/// Brackets have a level: `()` is 0, `[]` 1, `{}` 2, `<>` 3, and `Aa` to
/// `Zz` are 4 to 29. Writing a structure assigns the levels from 0 upward.
///
/// ```
/// use ff_structure::DotBracket;
///
/// let db = DotBracket::try_from('[').unwrap();
/// assert_eq!(db, DotBracket::open(1));
/// assert_eq!(db.open_level(), Some(1));
/// assert_eq!(db.close_level(), None);
/// assert_eq!(DotBracket::try_from('(').unwrap(), DotBracket::OPEN);
/// assert_eq!(DotBracket::try_from('&').unwrap(), DotBracket::STRAND_BREAK);
/// assert_eq!(char::from(DotBracket::STRAND_BREAK), '+');
/// assert_eq!(char::from(DotBracket::STRAND_WRAP), '^');
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct DotBracket(u8);

impl DotBracket {
    /// The number of bracket levels.
    pub const LEVELS: usize = 30;
    /// `.`
    pub const UNPAIRED: Self = DotBracket(b'.');
    /// `+`, also parsed from `&`: the strand before it ends.
    pub const STRAND_BREAK: Self = DotBracket(b'+');
    /// `^`: the strand before it is circular, it continues at its own first
    /// base.
    pub const STRAND_WRAP: Self = DotBracket(b'^');
    /// `(`
    pub const OPEN: Self = DotBracket::open(0);
    /// `)`
    pub const CLOSE: Self = DotBracket::close(0);

    /// The opening bracket of `level`. Panics if `level >= LEVELS`.
    pub const fn open(level: usize) -> Self {
        DotBracket(OPENING[level])
    }

    /// The closing bracket of `level`. Panics if `level >= LEVELS`.
    pub const fn close(level: usize) -> Self {
        DotBracket(CLOSING[level])
    }

    /// The level, if this is an opening bracket.
    pub fn open_level(self) -> Option<usize> {
        match self.0 {
            b'(' => Some(0),
            b'[' => Some(1),
            b'{' => Some(2),
            b'<' => Some(3),
            b'A'..=b'Z' => Some(4 + (self.0 - b'A') as usize),
            _ => None,
        }
    }

    /// The level, if this is a closing bracket.
    pub fn close_level(self) -> Option<usize> {
        match self.0 {
            b')' => Some(0),
            b']' => Some(1),
            b'}' => Some(2),
            b'>' => Some(3),
            b'a'..=b'z' => Some(4 + (self.0 - b'a') as usize),
            _ => None,
        }
    }
}

impl TryFrom<char> for DotBracket {
    type Error = StructureError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '+' | '&' => Ok(DotBracket::STRAND_BREAK),
            '^' => Ok(DotBracket::STRAND_WRAP),
            '.' | '(' | ')' | '[' | ']' | '{' | '}' | '<' | '>'
                | 'A'..='Z' | 'a'..='z' => Ok(DotBracket(c as u8)),
            _ => Err(StructureError::InvalidToken(
                format!("character '{}'", c), "dot-bracket".into(), 0)),
        }
    }
}

impl From<DotBracket> for char {
    fn from(db: DotBracket) -> Self {
        db.0 as char
    }
}

impl fmt::Debug for DotBracket {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", char::from(*self))
    }
}

/// A structure in dot-bracket notation, one byte per position. It is the
/// format to store many structures in, for example as keys of a hash map.
///
/// The inner `Vec` is public, and the positions can be changed without
/// checks:
///
/// ```
/// use ff_structure::{DotBracket, DotBracketVec};
///
/// let mut db = DotBracketVec::try_from("(...)").unwrap();
/// db[1] = DotBracket::OPEN;
/// db[3] = DotBracket::CLOSE;
/// assert_eq!(db.to_string(), "((.))");
/// ```
///
/// Equal structures can be written with different brackets, and then do
/// not compare equal. A conversion to a pair table and back gives one
/// representation for each structure:
///
/// ```
/// use ff_structure::{DotBracketVec, MultiPairTable};
///
/// let db = DotBracketVec::try_from("[(.)]").unwrap();
/// let mpt = MultiPairTable::try_from(&db).unwrap();
/// assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), "((.))");
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct DotBracketVec(pub Vec<DotBracket>);

impl Deref for DotBracketVec {
    type Target = [DotBracket];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for DotBracketVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/// Parse a `&str` into a `DotBracketVec`. An `InvalidToken` error carries 
/// the position of the character.
impl TryFrom<&str> for DotBracketVec {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut vec = Vec::with_capacity(s.len());
        for (i, c) in s.chars().enumerate() {
            match DotBracket::try_from(c) {
                Ok(db) => vec.push(db),
                Err(StructureError::InvalidToken(tok, src, _)) => {
                    return Err(StructureError::InvalidToken(tok, src, i));
                }
                Err(e) => return Err(e),
            }
        }
        Ok(DotBracketVec(vec))
    }
}

/// Levels are assigned greedily, left to right: each opening pair takes the
/// lowest level whose innermost open pair encloses it. All pairs of a band
/// are crossed by the same pairs, so each band gets one level; bands that do
/// not cross may share one.
impl TryFrom<&MultiPairTable> for DotBracketVec {
    type Error = StructureError;

    fn try_from(mpt: &MultiPairTable) -> Result<Self, Self::Error> {
        // Per level: the 3' ends of the currently open pairs, innermost on top.
        let mut open: Vec<Vec<usize>> = vec![Vec::new(); DotBracket::LEVELS];
        let mut level_of = vec![0usize; mpt.len()];
        let mut result = Vec::with_capacity(mpt.len());

        for (i, state) in mpt.iter().enumerate() {
            match state {
                MultiStruct::StrandBreak => result.push(DotBracket::STRAND_BREAK),
                MultiStruct::StrandWrap => result.push(DotBracket::STRAND_WRAP),
                MultiStruct::Unpaired => result.push(DotBracket::UNPAIRED),
                MultiStruct::Paired(j) if j > i => {
                    let k = open.iter()
                        .position(|s| s.last().is_none_or(|&q| j < q))
                        .ok_or(StructureError::TooManyPseudoknotLevels(i, j))?;
                    open[k].push(j);
                    level_of[i] = k;
                    result.push(DotBracket::open(k));
                }
                MultiStruct::Paired(j) => {
                    debug_assert!(j < i, "A base cannot pair with itself ({}).", i);
                    let k = level_of[j];
                    debug_assert_eq!(open[k].last(), Some(&i), "Pairs within a level must nest.");
                    open[k].pop();
                    result.push(DotBracket::close(k));
                }
            }
        }
        Ok(DotBracketVec(result))
    }
}

impl fmt::Display for DotBracketVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for db in &self.0 {
            write!(f, "{}", char::from(*db))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A table of `n` unpaired positions with the given pairs.
    fn table(n: usize, pairs: &[(usize, usize)]) -> MultiPairTable {
        let mut mpt = MultiPairTable::try_from(".".repeat(n).as_str()).unwrap();
        for &(i, j) in pairs {
            mpt.set(i, MultiStruct::Paired(j));
            mpt.set(j, MultiStruct::Paired(i));
        }
        mpt
    }

    #[test]
    fn test_write_back() {
        for (s, out) in [
            ("....", "...."),
            ("..(+.)..", "..(+.).."),
            ("((.+.)..)", "((.+.)..)"),
            ("([.+.](...)+...+)", "((.+.)(...)+...+)"),
            ("([{<)]}>", "([{<)]}>"),
            ("([{<>}])", "(((())))"),
            // Five crossing levels need a letter family.
            ("([{<A)]}>a", "([{<A)]}>a"),
            // Strand wraps: a circle, a circle bound by a linear strand, two
            // circles, a pseudoknot across a strand wrap.
            ("(((...)))^", "(((...)))^"),
            ("((((....^))))", "((((....^))))"),
            ("[[..]]^((..))^", "((..))^((..))^"),
            ("([.^.)]", "([.^.)]"),
        ] {
            let mpt = MultiPairTable::try_from(s).unwrap();
            assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), out, "{s}");
        }
    }

    #[test]
    fn test_levels_from_pairs() {
        // Each pair crosses all others: one level per pair.
        let mpt = table(6, &[(0, 3), (1, 4), (2, 5)]);
        assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), "([{)]}");
        // Pairs that do not cross share the lowest level.
        let mpt = table(8, &[(0, 7), (1, 2), (3, 5), (4, 6)]);
        assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), "(()([)])");
    }

    #[test]
    fn test_too_many_levels() {
        let crossing = |k: usize| table(2 * k, &(0..k).map(|i| (i, i + k)).collect::<Vec<_>>());
        let db = DotBracketVec::try_from(&crossing(DotBracket::LEVELS)).unwrap().to_string();
        assert_eq!(db, "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ)]}>abcdefghijklmnopqrstuvwxyz");
        let err = DotBracketVec::try_from(&crossing(DotBracket::LEVELS + 1)).unwrap_err();
        assert!(matches!(err, StructureError::TooManyPseudoknotLevels(30, 61)));
    }

    #[test]
    fn test_single_char_roundtrip() {
        for level in 0..DotBracket::LEVELS {
            let (open, close) = (DotBracket::open(level), DotBracket::close(level));
            assert_eq!(DotBracket::try_from(char::from(open)).unwrap(), open);
            assert_eq!(DotBracket::try_from(char::from(close)).unwrap(), close);
            assert_eq!(open.open_level(), Some(level));
            assert_eq!(close.close_level(), Some(level));
        }
        assert_eq!(DotBracket::try_from('.').unwrap(), DotBracket::UNPAIRED);
        assert_eq!(DotBracket::try_from('+').unwrap(), DotBracket::STRAND_BREAK);
        assert_eq!(DotBracket::try_from('&').unwrap(), DotBracket::STRAND_BREAK);
        assert_eq!(DotBracket::try_from('^').unwrap(), DotBracket::STRAND_WRAP);
    }

    #[test]
    fn test_vec_from_str() {
        let v = DotBracketVec::try_from("([.)]").unwrap();
        assert_eq!(v.0, vec![DotBracket::OPEN, DotBracket::open(1), DotBracket::UNPAIRED,
            DotBracket::CLOSE, DotBracket::close(1)]);
        let s = "([{<A.Bb.a>}])";
        assert_eq!(DotBracketVec::try_from(s).unwrap().to_string(), s);
    }

    #[test]
    fn test_one_byte() {
        assert_eq!(std::mem::size_of::<DotBracket>(), 1);
    }

    #[test]
    fn test_vec_invalid_token_position() {
        let err = DotBracketVec::try_from("(.#)").unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(_, ref src, 2) if src == "dot-bracket"));
    }
}
