use std::convert::TryFrom;
use std::ops::{Deref, DerefMut};
use crate::NAIDX;
use crate::StructureError;
use crate::DotBracket;
use crate::DotBracketVec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiPairTable(Vec<Vec<Option<(NAIDX, NAIDX)>>>);

impl MultiPairTable {
    /// Total number of nucleotides across all strands
    pub fn len(&self) -> usize {
        self.0.iter().map(|s| s.len()).sum()
    }

    /// True if there are no strands
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Number of strands
    pub fn num_strands(&self) -> usize {
        self.0.len()
    }

    // NOTE: we may provide a different utility function in the 
    // future that does usize <-> NAIDX conversions for readability.
    pub fn get_pair(&self, loc: (usize, usize)) -> &Option<(NAIDX, NAIDX)> {
        &self.0[loc.0][loc.1]
    }

}

impl Deref for MultiPairTable {
    type Target = [Vec<Option<(NAIDX, NAIDX)>>];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for MultiPairTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl TryFrom<&str> for MultiPairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut strand_idx = 0;
        let mut domain_idx = 0;
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (strand_idx, domain_idx)
        let mut pair_table: Vec<Vec<Option<(NAIDX, NAIDX)>>> = vec![vec![]];

        for (i, ch) in s.chars().enumerate() {
            match ch {
                '+' | '&' => {
                    if strand_idx == 0 && domain_idx == 0 {
                        return Err(StructureError::InvalidToken(
                                "strand break".into(),
                                "complex".into(),
                                0,
                        ));
                    }
                    // skip empty final strand.
                    if i < s.len()-1 {
                        pair_table.push(vec![]);
                    }
                    strand_idx += 1;
                    domain_idx = 0;
                }
                '(' => {
                    stack.push((strand_idx, domain_idx));
                    pair_table[strand_idx].push(None); // placeholder
                    domain_idx += 1;
                }
                ')' => {
                    let (si, di) = stack.pop()
                        .ok_or(StructureError::UnmatchedMultiClose((strand_idx, domain_idx)))?;
                    pair_table[si][di] = Some((strand_idx as NAIDX, domain_idx as NAIDX));
                    pair_table[strand_idx].push(Some((si as NAIDX, di as NAIDX)));
                    domain_idx += 1;
                }
                '.' => {
                    pair_table[strand_idx].push(None);
                    domain_idx += 1;
                }
                _ => {
                    return Err(StructureError::InvalidToken(
                        format!("character '{}'", ch), 
                        "complex".into(), 
                        i));
                }
            }
        }
        if let Some((si, di)) = stack.pop() {
            return Err(StructureError::UnmatchedMultiOpen((si, di)));
        }
        Ok(MultiPairTable(pair_table))
    }
}

impl TryFrom<&DotBracketVec> for MultiPairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        // Multi-stranded case
        let mut strand_idx = 0;
        let mut domain_idx = 0;
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (strand_idx, domain_idx)
        let mut pair_table: Vec<Vec<Option<(NAIDX, NAIDX)>>> = vec![vec![]];

        for dot in db.iter() {
            match dot {
                DotBracket::Break => {
                    if strand_idx == 0 && domain_idx == 0 {
                        return Err(StructureError::InvalidToken(
                                "strand break".into(),
                                "complex".into(),
                                0,
                        ));
                    }
                    pair_table.push(vec![]);
                    strand_idx += 1;
                    domain_idx = 0;
                }
                DotBracket::Open => {
                    stack.push((strand_idx, domain_idx));
                    pair_table[strand_idx].push(None); // placeholder
                    domain_idx += 1;
                }
                DotBracket::Close => {
                    let (si, di) = stack.pop()
                        .ok_or(StructureError::UnmatchedMultiClose(
                            (strand_idx, domain_idx)))?;
                    pair_table[si][di] = Some((strand_idx as NAIDX, domain_idx as NAIDX));
                    pair_table[strand_idx].push(Some((si as NAIDX, di as NAIDX)));
                    domain_idx += 1;
                }
                DotBracket::Unpaired => {
                    pair_table[strand_idx].push(None);
                    domain_idx += 1;
                }
            }
        }

        if let Some((si, di)) = stack.pop() {
            return Err(StructureError::UnmatchedMultiOpen((si, di)));
        }

        Ok(MultiPairTable(pair_table))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multi_pair_table() {
        let pt = MultiPairTable::try_from("((.+.))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((1, 2)));
        assert_eq!(*pt.get_pair((0, 1)), Some((1, 1)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((1, 0)), None);
        assert_eq!(*pt.get_pair((1, 1)), Some((0, 1)));
        assert_eq!(*pt.get_pair((1, 2)), Some((0, 0)));
    }

    #[test]
    fn test_multi_pair_table_hack() {
        let pt = MultiPairTable::try_from("((..))+").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((0, 5)));
        assert_eq!(*pt.get_pair((0, 1)), Some((0, 4)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((0, 3)), None);
        assert_eq!(*pt.get_pair((0, 4)), Some((0, 1)));
        assert_eq!(*pt.get_pair((0, 5)), Some((0, 0)));
    }
}


