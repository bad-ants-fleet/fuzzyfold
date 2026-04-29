//! PairList definition. 
//!
//! Compact integer-based representation of base pairs, can be used 
//! as alternative to PairTable and PairSet representations.
//!
//! A `Pair` is defined by two 16-bit indices (`NAIDX`) packed into a vector.
//! 
//! Example usage: when length of the exterior loop is not relevant,
//! e.g. in cotranscriptional simulations where structures with the same base pairs
//! are treated as the same structure independent of their exterior loop length. 
//!  

use std::fmt;
use std::ops::Deref;

use crate::DotBracket;
use crate::DotBracketVec;
use crate::PairTable;
use crate::Pair;
use crate::NAIDX;
use crate::StructureError;


/// All base pairs of a structre stored in a Vector, sorted by the index of the opening base
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct PairList {
    pairs: Vec<(NAIDX, NAIDX)>,
}

impl PairList {
    /// Create an empty pair list.
    pub fn new() -> Self {
        Self {
            pairs: Vec::new()
        }
    }

    /// Number of pairs contained in the list.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns true if there are no pairs.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Insert a new pair; returns true if it was newly inserted.
    pub fn insert(&mut self, pair: Pair) -> bool {
        let entry = (pair.i(), pair.j());
        if self.pairs.contains(&entry) {
            return false
        }else {
            self.pairs.push(entry);
            true
        }
    }

    /// Check if a pair exists in the list.
    pub fn contains(&self, pair: &Pair) -> bool {
        self.pairs.contains(&(pair.i(), pair.j()))
    }

}

impl Deref for PairList {
    type Target = Vec<(NAIDX, NAIDX)>;
    fn deref(&self) -> &Self::Target {
        &self.pairs
    }
}


impl From<&PairTable> for PairList {
    /// build PairList from PairTable 
    fn from(pt: &PairTable) -> Self { 
        let mut pairs = Vec::new();
        for (i, &j_opt) in pt.iter().enumerate() {
            let i = i as NAIDX;
            if let Some(j) = j_opt {
                if i < j {
                    pairs.push((i,j));
                }
            }
        }
        Self {
            pairs,
        }
    }
}


impl TryFrom<&DotBracketVec> for PairList {
    // build PairList from DotBrackVector
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut pairs = Vec::new();

        for (i, dot) in db.iter().enumerate() {
            match dot {
                DotBracket::Open => stack.push(i),
                DotBracket::Close => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    pairs.push((j as NAIDX, i as NAIDX));
                }
                _ => continue,
            }
        }
        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }

        pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        
        Ok(PairList {pairs})
    }
}



impl fmt::Display for PairList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for &(i,j) in &self.pairs {
            if !first {
                write!(f, ",")?;
            }
            // Only here we show 1-based values for readability.
            write!(f, "({},{})", i, j)?;
            first = false;
        }
        Ok(())
    }
}


#[cfg(test)]
mod tests {
    use super::*;


    #[test]
    fn test_pair_list_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::from(&pt);

        let expected = vec![Pair::new(0, 5), Pair::new(1, 4)];
        
        for p in &expected {
            assert!(pl.contains(p));
        }
        assert!(!pl.contains(&Pair::new(0, 4)));
    }

    #[test]
    fn test_pair_list_from_dot_bracket_vector() {
        let dbv = DotBracketVec::try_from("((..))").unwrap();
        let pl = PairList::try_from(&dbv).unwrap();

        assert_eq!(pl.len(), 2);
        assert!(pl.contains(&Pair::new(0, 5)));
        assert!(pl.contains(&Pair::new(1, 4)));
        assert!(!pl.contains(&Pair::new(0, 4)));
    }

    #[test]
    fn test_display() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::from(&pt);
        println!("PairList:{}", pl);
        let s = format!("{}", pl);
        assert!(s.contains("(0,5)"));
        assert!(s.contains("(1,4)"));
    }

    #[test]
    fn test_order() {
        let pt = PairTable::try_from("((..(...).))").unwrap();
        let pl = PairList::from(&pt);

        assert_eq!(pl[0], (0, 11));
        assert_eq!(pl[1], (1, 10));
    }
}