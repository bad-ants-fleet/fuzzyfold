//! Pair and PairList definitions. 
//!
//! They can be used as alternative to PairTable representations,
//! but beware that these implementations are 1-based. 
//! 
//! We currently do not povide the conversions from PairList to 
//! PairTable, etc., as sanity checks would be quite expensive,
//! and we don't want to produce invalid PairTables. 
//! Presumably, we will convert to DotBracketVec instead or work
//! with LoopIndex stuff.
//! 

use ff_structure::PairTable;
use ff_structure::NAIDX;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Pair { i: NAIDX, j: NAIDX }

impl Pair {
    pub fn new(i: NAIDX, j: NAIDX) -> Self {
        debug_assert!(i < j);
        Pair { i, j }
    }

    pub fn i(&self) -> NAIDX {
        self.i
    }

    pub fn j(&self) -> NAIDX {
        self.j
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairList {
    length: usize,
    pairs: Vec<Pair>,
}

impl PairList {
    pub fn pairs(&self) -> &Vec<Pair> {
        &self.pairs
    }
}

impl From<&PairTable> for PairList {
    fn from(pt: &PairTable) -> Self {
        let mut pairs = Vec::new();
        for (i, &j_opt) in pt.iter().enumerate() {
            let i = i as NAIDX;
            if let Some(j) = j_opt {
                if i > j {
                    continue;
                }
                pairs.push(Pair::new(i + 1, j + 1));
            }
        }
        PairList {
            length: pt.len(),
            pairs,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_list_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::from(&pt);

        assert_eq!(pl.length, 6);
        assert_eq!(pl.pairs, vec![Pair::new(1, 6), Pair::new(2, 5)]);
    }
}
