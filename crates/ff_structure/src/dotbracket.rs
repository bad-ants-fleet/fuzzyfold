use std::fmt;
use std::ops::Deref;
use std::ops::DerefMut;
use std::convert::TryFrom;

use crate::PairTable;
use crate::MultiPairTable;
use crate::StructureError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DotBracket {
    Unpaired, // '.'
    Open,     // '('
    Close,    // ')'
    Break,    // '+' or '&'
}

impl TryFrom<char> for DotBracket {
    type Error = StructureError;

    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '.' => Ok(DotBracket::Unpaired),
            '(' => Ok(DotBracket::Open),
            ')' => Ok(DotBracket::Close),
            '+' | '&' => Ok(DotBracket::Break),
            _ => Err(StructureError::InvalidToken(c.to_string(), "dot-bracket".into(), 0)),
        }
    }
}

impl From<DotBracket> for char {
    fn from(db: DotBracket) -> Self {
        match db {
            DotBracket::Open => '(',
            DotBracket::Close => ')',
            DotBracket::Unpaired => '.',
            DotBracket::Break => '+',
        }
    }
}

/// DotBracketVec is a compact representation of secondary structure. Note that
/// the field is public, to allow unsafe modifications. Thus, DotBracketVecs can
/// be malformed and should be converted using the TryFrom trait.
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

impl From<&PairTable> for DotBracketVec {
    fn from(pt: &PairTable) -> Self {
        let mut result: Vec<DotBracket> = Vec::with_capacity(pt.len());
        for (i, &j_opt) in pt.iter().enumerate() {
            match j_opt {
                None => result.push(DotBracket::Unpaired),
                Some(j) if (j as usize) > i => result.push(DotBracket::Open),
                Some(j) if (j as usize) < i => result.push(DotBracket::Close),
                Some(j) if (j as usize) == i => {
                    unreachable!("PairTable construction prevents self-pairing! ({})", i);
                }
                _ => unreachable!(),
            }
        }
        DotBracketVec(result)
    }
}


impl From<&MultiPairTable> for DotBracketVec {

    fn from(pt: &MultiPairTable) -> Self {
        let mut result: Vec<DotBracket> = Vec::with_capacity(pt.len() + pt.num_strands());

        for (si, strand) in pt.iter().enumerate() {
            for (di, &pair) in strand.iter().enumerate() {
                match pair {
                    None => result.push(DotBracket::Unpaired),
                    Some((sj, dj)) => {
                        let sj = sj as usize;
                        let dj = dj as usize;
                        if (sj, dj) > (si, di) {
                            result.push(DotBracket::Open);
                        } else if (sj, dj) < (si, di) {
                            result.push(DotBracket::Close);
                        } else {
                            panic!("Invalid self-pairing at strand {si}, domain {di}");
                        }
                    }
                }
            }
            // NOTE:: pushes strand break at the end of the DotBracketVec,
            // intentionally!
            result.push(DotBracket::Break);
        }

        DotBracketVec(result)
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

    #[test]
    fn test_dot_bracket_from_char() {
        assert_eq!(DotBracket::try_from('.').unwrap(), DotBracket::Unpaired);
        assert_eq!(DotBracket::try_from('(').unwrap(), DotBracket::Open);
        assert_eq!(DotBracket::try_from(')').unwrap(), DotBracket::Close);
    }

    #[test]
    fn test_char_from_dot_bracket() {
        assert_eq!(char::from(DotBracket::Unpaired), '.');
        assert_eq!(char::from(DotBracket::Open), '(');
        assert_eq!(char::from(DotBracket::Close), ')');
    }

    #[test]
    fn test_dot_bracket_from_invalid_char() {
        let res = DotBracket::try_from('x');
        assert!(matches!(res, Err(StructureError::InvalidToken(_, src, _)) if src == "dot-bracket"));
    }

    #[test]
    fn test_dot_bracket_vec_from_str() {
        let dbv = DotBracketVec::try_from("(.).").unwrap();
        assert_eq!(format!("{}", dbv), "(.).");
        assert_eq!(dbv.len(), 4);
        assert_eq!(dbv[0], DotBracket::Open);
        assert_eq!(dbv[1], DotBracket::Unpaired);
        assert_eq!(dbv[2], DotBracket::Close);
        assert_eq!(dbv[3], DotBracket::Unpaired);
    }

    #[test]
    fn test_dot_bracket_vec_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let dbv = DotBracketVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((..))");
    }

    #[test]
    fn test_dot_bracket_vec_from_multi_pair_table_hack() {
        let pt = MultiPairTable::try_from("((..))+").unwrap();
        let dbv = DotBracketVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((..))+");
    }

    #[test]
    fn test_dot_bracket_vec_from_multi_pair_table() {
        let pt = MultiPairTable::try_from("((..)+)").unwrap();
        let dbv = DotBracketVec::from(&pt);
        assert_eq!(format!("{}", dbv), "((..)+)+");
    }

}

