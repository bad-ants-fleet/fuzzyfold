//! Pair and PairSet definitions. 
//!
//! Compact integer-based representation of base pairs, can 
//! be used as alternative to PairTable representations.
//!
//! A `Pair` is defined by two 16-bit indices (`NAIDX`) packed into a
//! 32-bit integer key for efficient set and map storage.
//!
//! Indices are **0-based** throughout, consistent with `PairTable`,
//! `Array2` matrices, and Nussinov dynamic programming indices.
//!
//! We currently do not povide the conversions from PairSet to 
//! PairTable, etc., as sanity checks would be quite expensive,
//! and we don't want to produce invalid PairTables. 
//! Presumably, we will convert to DotBracketVec instead or work
//! with LoopIndex stuff.
//! 

use std::fmt;
use nohash_hasher::IntSet;
use ff_structure::PairTable;
use ff_structure::NAIDX;

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
        //debug_assert!(i < (1 << 16));
        //debug_assert!(j < (1 << 16));
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairSet {
    length: usize,
    pairs: IntSet<P1KEY>,
}

impl PairSet {
    /// Create an empty pair set for a given sequence length.
    pub fn new(length: usize) -> Self {
        Self {
            length,
            pairs: IntSet::default(),
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
        debug_assert!((pair.j() as usize) < self.length);
        self.pairs.insert(pair.key())
    }

    /// Check if a pair exists in the set.
    pub fn contains(&self, pair: &Pair) -> bool {
        self.pairs.contains(&pair.key())
    }

    /// Iterator over all pairs in arbitrary order.
    pub fn iter(&self) -> impl Iterator<Item = Pair> + '_ {
        self.pairs.iter().map(|&k| Pair::from_key(k))
    }

    /// Return all pairs as a Vec (for deterministic inspection).
    pub fn to_vec(&self) -> Vec<Pair> {
        let mut v: Vec<_> = self.iter().collect();
        v.sort_unstable_by_key(|p| (p.i(), p.j()));
        v
    }

    /// Underlying sequence length (from the originating `PairTable`).
    pub fn length(&self) -> usize {
        self.length
    }
}

impl From<&PairTable> for PairSet {
    fn from(pt: &PairTable) -> Self {
        let mut pairs = IntSet::default();
        for (i, &j_opt) in pt.iter().enumerate() {
            let i = i as NAIDX;
            if let Some(j) = j_opt {
                if i < j {
                    pairs.insert(Pair::new(i, j).key());
                }
            }
        }
        Self {
            length: pt.len(),
            pairs,
        }
    }
}

impl fmt::Display for PairSet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut first = true;
        for pair in self.to_vec() {
            if !first {
                write!(f, ",")?;
            }
            // Only here we show 1-based values for readability.
            write!(f, "({},{})", pair.i(), pair.j())?;
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
        let pl = PairSet::from(&pt);

        let expected = vec![Pair::new(0, 5), Pair::new(1, 4)];
        assert_eq!(pl.length(), 6);
        assert_eq!(pl.to_vec(), expected);

        for p in &expected {
            assert!(pl.contains(p));
        }
        assert!(!pl.contains(&Pair::new(0, 4)));
    }

    #[test]
    fn test_display() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairSet::from(&pt);
        let s = format!("{}", pl);
        assert!(s.contains("(0,5)"));
        assert!(s.contains("(1,4)"));
    }
}

