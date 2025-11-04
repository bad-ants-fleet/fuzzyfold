//! PairTable construction and helper traits.

use std::ops::{Deref, DerefMut};
use std::convert::TryFrom;
use crate::NAIDX;
use crate::StructureError;
use crate::{DotBracket, DotBracketVec};

/// As of v0.1.3 the PairTabel field is private. A pair-table should
/// be constructed by From or TryFrom traits, but then be save to use.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairTable(Vec<Option<NAIDX>>);

impl PairTable {
    /// Check if the substructure from `i..j` is well-formed:
    /// - All pairings are internal to the interval
    pub fn is_well_formed(&self, i: usize, j: usize) -> bool {
        assert!(j <= self.len(), "Invalid interval: j must be <= length");

        for k in i..j {
            if let Some(l) = self[k] {
                let ul = l as usize;
                if ul < i || ul >= j {
                    return false; // points outside
                }
            }
        }
        true
    }
}

impl Deref for PairTable {
    type Target = [Option<NAIDX>];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for PairTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl TryFrom<&str> for PairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut table = vec![None; s.len()];

        for (i, c) in s.chars().enumerate() {
            match c {
                '(' => stack.push(i),
                ')' => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    table[i] = Some(j as NAIDX);
                    table[j] = Some(i as NAIDX);
                }
                '.' => (),
                _ => return Err(StructureError::InvalidToken(format!("character '{}'", c), "structure".to_string(), i)),
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }
        Ok(PairTable(table))
    }
}

impl TryFrom<&DotBracketVec> for PairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut table = vec![None; db.len()];

        for (i, dot) in db.iter().enumerate() {
            match dot {
                DotBracket::Open => stack.push(i),
                DotBracket::Close => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    table[i] = Some(j as NAIDX);
                    table[j] = Some(i as NAIDX);
                }
                DotBracket::Unpaired => {}
                DotBracket::Break => unreachable!("unexpected Break in single-stranded case"),
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        Ok(PairTable(table))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(pt[0], Some(5));
        assert_eq!(pt[1], Some(4));
        assert_eq!(pt[2], None);
        assert_eq!(pt[3], None);
        assert_eq!(pt[4], Some(1));
        assert_eq!(pt[5], Some(0));
    }

    #[test]
    fn test_unmatched_open() {
        let err = PairTable::try_from("(()").unwrap_err();
        assert_eq!(format!("{}", err), "Unmatched '(' at position 0");
    }

    #[test]
    fn test_unmatched_close() {
        let err = PairTable::try_from("())").unwrap_err();
        assert_eq!(format!("{}", err), "Unmatched ')' at position 2");
    }

    #[test]
    fn test_invalid_token() {
        let err = PairTable::try_from("(x)").unwrap_err();
        assert_eq!(format!("{}", err), "Invalid character 'x' in structure at position 1");
    }

    #[test]
    fn test_well_formed_empty_interval() {
        let pt= PairTable::try_from("...").unwrap();
        assert!(pt.is_well_formed(0, 0)); 
        assert!(pt.is_well_formed(0, 1)); 
        assert!(pt.is_well_formed(0, 2)); 
        assert!(pt.is_well_formed(0, 3)); 
        assert!(pt.is_well_formed(1, 3)); 
        assert!(pt.is_well_formed(2, 3)); 
        assert!(pt.is_well_formed(3, 3)); 
    }

    #[test]
    fn test_well_formed_pairings_within_interval() {
        let pt = PairTable::try_from(".(.).").unwrap();
        assert!(pt.is_well_formed(0, 5)); // Full interval -- 0-based
        assert!(pt.is_well_formed(0, 4)); 
        assert!(pt.is_well_formed(1, 5));
        assert!(pt.is_well_formed(1, 4));
        assert!(pt.is_well_formed(1, 4));
        assert!(pt.is_well_formed(2, 3));
        assert!(!pt.is_well_formed(0, 3)); 
        assert!(!pt.is_well_formed(1, 3)); 
        assert!(!pt.is_well_formed(2, 4)); 
    }

    #[test]
    #[should_panic(expected = "Invalid interval: j must be <= length")]
    fn test_well_formed_out_of_bounds_assert() {
        let pt = PairTable::try_from("..").unwrap();
        pt.is_well_formed(0, 3); // j = pt.len(), should panic
    }
}



