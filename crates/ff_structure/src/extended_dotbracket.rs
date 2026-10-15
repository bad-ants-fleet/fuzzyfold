//! Extended-dot-bracket notation.
//!
//! Keeps the canonical dot-bracket implementation separate from the
//! extended-dot-bracket. Especially for pseudoknots, we need more characters
//! for nested structures. This implementation uses two separate types:
//! - `ExtendedDotBracket`  — a single token
//! - `ExtendedDotBracketVec` — a sequence of tokens (mirrors `DotBracketVec`)

use std::fmt;
use std::ops::{Deref, DerefMut};
use std::convert::TryFrom;

use crate::DotBracket;
use crate::DotBracketVec;
use crate::PairTable;
use crate::StructureError;

// ---------------------------------------------------------------------------
// BracketKind
// ---------------------------------------------------------------------------

/// The bracket type, ordered by conventional pseudoknot nesting level.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BracketKind {
    Round,   // ( )
    Square,  // [ ]
    Curly,   // { }
    Angle,   // < >
    /// Letter families, `0 => A a` through `25 => Z z`.
    Letter(u8),
}

/// Number of letter families, `Aa`..`Zz`.
const N_LETTERS: u8 = 26;

/// `([{<` followed by `Aa`..`Zz`, the conventional pseudoknot-nesting order.
static ALL: [BracketKind; 4 + N_LETTERS as usize] = {
    let mut all = [BracketKind::Round; 4 + N_LETTERS as usize];
    all[1] = BracketKind::Square;
    all[2] = BracketKind::Curly;
    all[3] = BracketKind::Angle;
    let mut k = 0;
    while k < N_LETTERS {
        all[4 + k as usize] = BracketKind::Letter(k);
        k += 1;
    }
    all
};

impl BracketKind {
    pub fn open_char(self) -> char {
        match self {
            BracketKind::Round     => '(',
            BracketKind::Square    => '[',
            BracketKind::Curly     => '{',
            BracketKind::Angle     => '<',
            BracketKind::Letter(k) => (b'A' + k) as char,
        }
    }

    pub fn close_char(self) -> char {
        match self {
            BracketKind::Round     => ')',
            BracketKind::Square    => ']',
            BracketKind::Curly     => '}',
            BracketKind::Angle     => '>',
            BracketKind::Letter(k) => (b'a' + k) as char,
        }
    }

    /// All kinds in their conventional pseudoknot-nesting order.
    pub fn all() -> &'static [BracketKind] {
        &ALL
    }
}

// ---------------------------------------------------------------------------
// ExtendedDotBracket — single token
// ---------------------------------------------------------------------------

/// Represents a single character in extended dot-bracket notation,
/// covering simple and pseudoknot-encoding bracket types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExtendedDotBracket {
    Unpaired,           // '.'
    Break,              // '+' or '&'
    Open(BracketKind),
    Close(BracketKind),
}

impl TryFrom<char> for ExtendedDotBracket {
    type Error = StructureError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        use BracketKind::*;
        use ExtendedDotBracket::{Break, Close, Open, Unpaired};
        Ok(match c {
            '.' => Unpaired,
            '+' | '&' => Break,
            '(' => Open(Round),
            ')' => Close(Round),
            '[' => Open(Square),
            ']' => Close(Square),
            '{' => Open(Curly),
            '}' => Close(Curly),
            '<' => Open(Angle),
            '>' => Close(Angle),
            'A'..='Z' => Open(Letter(c as u8 - b'A')),
            'a'..='z' => Close(Letter(c as u8 - b'a')),
            _ => {
                return Err(StructureError::InvalidToken(
                    format!("character '{}'", c),
                    "extended dot-bracket".to_string(),
                    0,
                ));
            }
        })
    }
}

impl From<ExtendedDotBracket> for char {
    fn from(edb: ExtendedDotBracket) -> Self {
        match edb {
            ExtendedDotBracket::Unpaired    => '.',
            ExtendedDotBracket::Break       => '+',
            ExtendedDotBracket::Open(kind)  => kind.open_char(),
            ExtendedDotBracket::Close(kind) => kind.close_char(),
        }
    }
}

/// Lossless downcast: only succeeds when only Round brackets are present.
impl TryFrom<ExtendedDotBracket> for DotBracket {
    type Error = StructureError;

    fn try_from(edb: ExtendedDotBracket) -> Result<Self, Self::Error> {
        match edb {
            ExtendedDotBracket::Unpaired                  => Ok(DotBracket::Unpaired),
            ExtendedDotBracket::Break                     => Ok(DotBracket::Break),
            ExtendedDotBracket::Open(BracketKind::Round)  => Ok(DotBracket::Open),
            ExtendedDotBracket::Close(BracketKind::Round) => Ok(DotBracket::Close),
            other => Err(StructureError::InvalidToken(
                format!("character '{}'", char::from(other)),
                "simple dot-bracket".to_string(),
                0,
            )),
        }
    }
}

// ---------------------------------------------------------------------------
// ExtendedDotBracketVec — sequence of tokens (mirrors DotBracketVec)
// ---------------------------------------------------------------------------

/// A sequence of `ExtendedDotBracket` tokens. The inner `Vec` is public to
/// allow unsafe modifications, mirroring `DotBracketVec`.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ExtendedDotBracketVec(pub Vec<ExtendedDotBracket>);

impl Deref for ExtendedDotBracketVec {
    type Target = [ExtendedDotBracket];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for ExtendedDotBracketVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

/// Parse a `&str` into an `ExtendedDotBracketVec`, propagating position
/// information into any `InvalidToken` error — mirrors `DotBracketVec`.
impl TryFrom<&str> for ExtendedDotBracketVec {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut vec = Vec::with_capacity(s.len());
        for (i, c) in s.chars().enumerate() {
            match ExtendedDotBracket::try_from(c) {
                Ok(edb) => vec.push(edb),
                Err(StructureError::InvalidToken(tok, src, _)) => {
                    return Err(StructureError::InvalidToken(tok, src, i));
                }
                Err(e) => return Err(e),
            }
        }
        Ok(ExtendedDotBracketVec(vec))
    }
}

/// Lossless upcast: every `DotBracketVec` is a valid `ExtendedDotBracketVec`.
impl From<&DotBracketVec> for ExtendedDotBracketVec {
    fn from(dbv: &DotBracketVec) -> Self {
        let vec = dbv.iter().map(|&db| match db {
            DotBracket::Unpaired => ExtendedDotBracket::Unpaired,
            DotBracket::Break    => ExtendedDotBracket::Break,
            DotBracket::Open     => ExtendedDotBracket::Open(BracketKind::Round),
            DotBracket::Close    => ExtendedDotBracket::Close(BracketKind::Round),
        }).collect();
        ExtendedDotBracketVec(vec)
    }
}

/// Lossless downcast: only succeeds when no non-Round brackets are present.
impl TryFrom<&ExtendedDotBracketVec> for DotBracketVec {
    type Error = StructureError;

    fn try_from(edbv: &ExtendedDotBracketVec) -> Result<Self, Self::Error> {
        let vec = edbv.iter()
            .map(|&edb| DotBracket::try_from(edb))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(DotBracketVec(vec))
    }
}

/// Build from a `PairTable` — delegates through `DotBracketVec` since a plain
/// `PairTable` has no pseudoknot information.
impl From<&PairTable> for ExtendedDotBracketVec {
    fn from(pt: &PairTable) -> Self {
        ExtendedDotBracketVec::from(&DotBracketVec::from(pt))
    }
}

impl fmt::Display for ExtendedDotBracketVec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for edb in &self.0 {
            write!(f, "{}", char::from(*edb))?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_char_roundtrip() {
        for &kind in BracketKind::all() {
            let open  = ExtendedDotBracket::Open(kind);
            let close = ExtendedDotBracket::Close(kind);
            assert_eq!(ExtendedDotBracket::try_from(char::from(open)).unwrap(),  open);
            assert_eq!(ExtendedDotBracket::try_from(char::from(close)).unwrap(), close);
        }
        assert_eq!(ExtendedDotBracket::try_from('.').unwrap(), ExtendedDotBracket::Unpaired);
        assert_eq!(ExtendedDotBracket::try_from('+').unwrap(), ExtendedDotBracket::Break);
        assert_eq!(ExtendedDotBracket::try_from('&').unwrap(), ExtendedDotBracket::Break);
    }

    #[test]
    fn test_invalid_char() {
        let err = ExtendedDotBracket::try_from('?').unwrap_err();
        assert!(matches!(
            err,
            StructureError::InvalidToken(_, ref src, _) if src == "extended dot-bracket"
        ));
    }

    #[test]
    fn test_vec_from_str_simple() {
        let v = ExtendedDotBracketVec::try_from("(.)").unwrap();
        assert_eq!(v.len(), 3);
        assert_eq!(v[0], ExtendedDotBracket::Open(BracketKind::Round));
        assert_eq!(v[1], ExtendedDotBracket::Unpaired);
        assert_eq!(v[2], ExtendedDotBracket::Close(BracketKind::Round));
        assert_eq!(format!("{}", v), "(.)");
    }

    #[test]
    fn test_vec_from_str_pseudoknot() {
        let v = ExtendedDotBracketVec::try_from("([)]").unwrap();
        assert_eq!(v[0], ExtendedDotBracket::Open(BracketKind::Round));
        assert_eq!(v[1], ExtendedDotBracket::Open(BracketKind::Square));
        assert_eq!(v[2], ExtendedDotBracket::Close(BracketKind::Round));
        assert_eq!(v[3], ExtendedDotBracket::Close(BracketKind::Square));
        assert_eq!(format!("{}", v), "([)]");
    }

    /// Every kind, opened in order and closed in reverse, so the string tracks
    /// `all()` instead of pinning a fixed alphabet.
    #[test]
    fn test_vec_from_str_all_kinds() {
        let kinds = BracketKind::all();
        let s: String = kinds
            .iter()
            .map(|k| k.open_char())
            .chain(kinds.iter().rev().map(|k| k.close_char()))
            .collect();
        let v = ExtendedDotBracketVec::try_from(s.as_str()).unwrap();
        assert_eq!(format!("{}", v), s);
    }

    #[test]
    fn test_letter_families_reach_z() {
        assert_eq!(BracketKind::all().len(), 30);
        assert_eq!(
            ExtendedDotBracket::try_from('Z').unwrap(),
            ExtendedDotBracket::Open(BracketKind::Letter(25))
        );
        assert_eq!(
            ExtendedDotBracket::try_from('z').unwrap(),
            ExtendedDotBracket::Close(BracketKind::Letter(25))
        );
        let v = ExtendedDotBracketVec::try_from("(Zz)").unwrap();
        assert_eq!(format!("{}", v), "(Zz)");
    }

    #[test]
    fn test_vec_invalid_token_position() {
        let err = ExtendedDotBracketVec::try_from("(.?)").unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(_, _, 2)));
    }

    #[test]
    fn test_upcast_from_dot_bracket_vec() {
        let dbv = DotBracketVec::try_from("((..))").unwrap();
        let edbv = ExtendedDotBracketVec::from(&dbv);
        assert_eq!(format!("{}", edbv), "((..))");
    }

    #[test]
    fn test_downcast_to_dot_bracket_vec_ok() {
        let edbv = ExtendedDotBracketVec::try_from("((..))").unwrap();
        let dbv = DotBracketVec::try_from(&edbv).unwrap();
        assert_eq!(format!("{}", dbv), "((..))");
    }

    #[test]
    fn test_downcast_to_dot_bracket_vec_fails_on_pseudoknot() {
        let edbv = ExtendedDotBracketVec::try_from("([)]").unwrap();
        assert!(DotBracketVec::try_from(&edbv).is_err());
    }

    #[test]
    fn test_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let edbv = ExtendedDotBracketVec::from(&pt);
        assert_eq!(format!("{}", edbv), "((..))");
    }
}