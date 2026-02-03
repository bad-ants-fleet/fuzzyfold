//! Pair and PairList definitions. 
//!
//! Compact integer-based representation of base pairs, can 
//! be used as alternative to PairTable representations.
//!
//! A `Pair` is defined by two 16-bit indices (`NAIDX`) packed into a
//! 32-bit integer key (`P1KEY`) for efficient set and map storage.
//!
//! We currently do not povide the conversions from PairList to 
//! PairTable, mainly because at this stage it is not clear if
//! PairSet may be used in the future to include pseudoknots. 
//! 

use std::fmt;
use std::ops::Deref;

use crate::PairTable;
use crate::NAIDX;
use crate::P1KEY;


/// A base pair (i, j) with i < j.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Pair {
    i: NAIDX,
    j: NAIDX,
}

impl Pair {
    /// Create a new pair (i, j). Panics in debug if i >= j.
    pub fn new(i: NAIDX, j: NAIDX) -> Self {
        debug_assert!(i < j);
        debug_assert!(j < NAIDX::MAX);
        Pair { i, j }
    }

    /// Return the 5'-side index.
    pub fn i(&self) -> NAIDX {
        self.i
    }

    /// Return the 3'-side index.
    pub fn j(&self) -> NAIDX {
        self.j
    }

    /// Compact 32-bit key encoding both indices.
    pub fn key(&self) -> P1KEY {
        ((self.i as P1KEY) << 16) | (self.j as P1KEY)
    }

    /// Decode a key back into a `Pair`.
    pub fn from_key(key: P1KEY) -> Self {
        let i = (key >> 16) as NAIDX;
        let j = (key & 0xFFFF) as NAIDX;
        debug_assert!(i < j);
        Pair { i, j }
    }
}

/// A collection of base pairs represented as compact integer keys.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Default)]
pub struct PairList {
    pairs: Vec<(NAIDX, NAIDX)>,
}

impl PairList {
    /// Create an empty pair set for a given sequence length.
    pub fn new() -> Self {
        Self {
            pairs: Vec::new()
        }
    }

    /// Number of pairs contained in the set.
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

    /// Check if a pair exists in the set.
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
    fn test_pair_key_roundtrip() {
        let p = Pair::new(1, 42);
        let k = p.key();
        let q = Pair::from_key(k);
        assert_eq!(p, q);
    }

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
    fn test_display() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::from(&pt);
        println!("PairList:{}", pl);
        let s = format!("{}", pl);
        assert!(s.contains("(0,5)"));
        assert!(s.contains("(1,4)"));
    }
}

