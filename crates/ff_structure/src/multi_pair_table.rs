use std::convert::TryFrom;
use std::ops::{Deref, DerefMut};
use crate::NAIDX;
use crate::StructureError;
use crate::DotBracket;
use crate::DotBracketVec;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MultiStruct {
    Paired(NAIDX),
    Unpaired,
    StrandBreak,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiPairTable(Vec<MultiStruct>);

impl Deref for MultiPairTable {
    type Target = [MultiStruct];
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
        let mut stack = Vec::new();
        let mut table = vec![MultiStruct::Unpaired; s.len()];

        for (i, c) in s.chars().enumerate() {
            match c {
                '(' => stack.push(i),
                ')' => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    table[i] = MultiStruct::Paired(j as NAIDX);
                    table[j] = MultiStruct::Paired(i as NAIDX);
                }
                '.' => (),
                '+' | '&' =>  table[i] = MultiStruct::StrandBreak,
                _ => return Err(StructureError::InvalidToken(format!("character '{}'", c), "structure".to_string(), i)),
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }
        Ok(MultiPairTable(table))
    }
}

impl TryFrom<&DotBracketVec> for MultiPairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        let mut stack = Vec::new();
        let mut table = vec![MultiStruct::Unpaired; db.len()];

        for (i, c) in db.iter().enumerate() {
            match c {
                DotBracket::Open => stack.push(i),
                DotBracket::Close => {
                    let j = stack.pop().ok_or(StructureError::UnmatchedClose(i))?;
                    table[i] = MultiStruct::Paired(j as NAIDX);
                    table[j] = MultiStruct::Paired(i as NAIDX);
                }
                DotBracket::Unpaired => (),
                DotBracket::Break =>  table[i] = MultiStruct::StrandBreak,
            }
        }

        if let Some(i) = stack.pop() {
            return Err(StructureError::UnmatchedOpen(i));
        }
        Ok(MultiPairTable(table))
    }
}

impl From<&StrandPairTable> for MultiPairTable {
    fn from(spt: &StrandPairTable) -> Self {
        // Step 1: compute prefix offsets for each strand
        // offset[strand] = global index where that strand starts
        let mut offsets = Vec::with_capacity(spt.num_strands());
        let mut acc = 0usize;
        for strand in spt.iter() {
            offsets.push(acc);
            acc += strand.len() + 1; // +1 for StrandBreak
        }

        // Step 2: build flattened table
        let mut flat = Vec::with_capacity(acc);

        for strand in spt.iter() {
            for entry in strand.iter() {
                match entry {
                    Some((s, p)) => {
                        let partner_global = offsets[*s as usize] + (*p as usize);
                        flat.push(MultiStruct::Paired(partner_global as NAIDX));
                    }
                    None => flat.push(MultiStruct::Unpaired),
                }
            }

            // Insert break after each strand except the last
            flat.push(MultiStruct::StrandBreak);
        }

        // Remove final trailing StrandBreak (canonical form)
        if let Some(MultiStruct::StrandBreak) = flat.last() {
            flat.pop();
        }

        MultiPairTable(flat)
    }
}


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StrandPairTable(Vec<Vec<Option<(NAIDX, NAIDX)>>>);

impl StrandPairTable {
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

impl Deref for StrandPairTable {
    type Target = [Vec<Option<(NAIDX, NAIDX)>>];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for StrandPairTable {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl TryFrom<&str> for StrandPairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut strand_idx = 0;
        let mut domain_idx = 0;
        let mut stack: Vec<(usize, usize)> = Vec::new(); // (strand_idx, domain_idx)
        let mut pair_table: Vec<Vec<Option<(NAIDX, NAIDX)>>> = vec![vec![]];

        for (i, ch) in s.chars().enumerate() {
            match ch {
                '+' | '&' => {
                    if domain_idx == 0 {
                        return Err(StructureError::InvalidToken(
                                "strand break".into(),
                                "complex".into(),
                                strand_idx,
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
        Ok(StrandPairTable(pair_table))
    }
}

impl TryFrom<&DotBracketVec> for StrandPairTable {
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

        Ok(StrandPairTable(pair_table))
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_pair_table() {
        let pt = StrandPairTable::try_from("((.+.))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((1, 2)));
        assert_eq!(*pt.get_pair((0, 1)), Some((1, 1)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((1, 0)), None);
        assert_eq!(*pt.get_pair((1, 1)), Some((0, 1)));
        assert_eq!(*pt.get_pair((1, 2)), Some((0, 0)));
    }

    #[test]
    fn test_strand_pair_table_hack() {
        let pt = StrandPairTable::try_from("((..))+").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(*pt.get_pair((0, 0)), Some((0, 5)));
        assert_eq!(*pt.get_pair((0, 1)), Some((0, 4)));
        assert_eq!(*pt.get_pair((0, 2)), None);
        assert_eq!(*pt.get_pair((0, 3)), None);
        assert_eq!(*pt.get_pair((0, 4)), Some((0, 1)));
        assert_eq!(*pt.get_pair((0, 5)), Some((0, 0)));
    }

    #[test]
    fn test_multi_pair_table_empty() {
        let inp = "....";
        let spt = StrandPairTable::try_from(inp).unwrap();
        let mpt = MultiPairTable::from(&spt);
        let dbr = DotBracketVec::from(&mpt);
        println!("{}", dbr);
        assert_eq!(inp, dbr.to_string());
    }

    #[test]
    fn test_multi_pair_table_single() {
        let inp = "..(+.)..";
        let spt = StrandPairTable::try_from(inp).unwrap();
        let mpt = MultiPairTable::from(&spt);
        let dbr = DotBracketVec::from(&mpt);
        println!("{}", dbr);
        assert_eq!(inp, dbr.to_string());
    }

    #[test]
    fn test_multi_pair_table_multi() {
        let inp = "((.+.)(...)+...+)";
        let spt = StrandPairTable::try_from(inp).unwrap();
        let mpt = MultiPairTable::from(&spt);
        let dbr = DotBracketVec::from(&mpt);
        println!("{:?}", mpt);
        println!("{}", dbr);
        assert_eq!(inp, dbr.to_string());
    }

    #[test]
    fn test_multi_pair_table_invalid_01() {
        let inp = "((.++.)..)";
        let err = StrandPairTable::try_from(inp).unwrap_err();
        match err {
            StructureError::InvalidToken(what, ctx, pos) => {
                assert_eq!(what, "strand break");
                assert_eq!(ctx, "complex");
                assert_eq!(pos, 1);
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }


    #[test]
    fn test_multi_pair_table_invalid_02() {
        let inp = "+..)";
        let err = StrandPairTable::try_from(inp).unwrap_err();

        match err {
            StructureError::InvalidToken(what, ctx, pos) => {
                assert_eq!(what, "strand break");
                assert_eq!(ctx, "complex");
                assert_eq!(pos, 0);
            }
            other => panic!("unexpected error: {:?}", other),
        }
    }

}


