
use std::convert::TryFrom;
use ff_structure::PairTable;
use ff_structure::StructureError;
use ff_structure::NAIDX;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Pair(pub NAIDX, pub NAIDX);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairList {
    pub length: usize,
    pub pairs: Vec<Pair>,
}

impl TryFrom<&PairTable> for PairList {
    type Error = StructureError;

    fn try_from(pt: &PairTable) -> Result<Self, Self::Error> {
        let mut pairs = Vec::new();
        for (i, &j_opt) in pt.iter().enumerate() {
            if let Some(j) = j_opt {
                let j = j as usize;
                if j == i {
                    return Err(StructureError::InvalidToken("self-pairing".to_string(), "pair table".to_string(), i));
                }
                if j >= pt.len() {
                    return Err(StructureError::UnmatchedOpen(i));
                }
                match pt.get(j) {
                    Some(&Some(k)) if k as usize == i => {
                        if i < j {
                            pairs.push(Pair((i + 1) as NAIDX, (j + 1) as NAIDX));
                        }
                    }
                    _ => return Err(StructureError::UnmatchedOpen(i)),
                }
            }
        }
        Ok(PairList {
            length: pt.len(),
            pairs,
        })
    }
}

impl TryFrom<&PairList> for PairTable {
    type Error = StructureError;

    fn try_from(pl: &PairList) -> Result<Self, Self::Error> {
        let mut table = vec![None; pl.length];

        for &Pair(i1, j1) in &pl.pairs {
            let i1 = i1 as usize;
            let j1 = j1 as usize;
            // Convert 1-based to 0-based
            if i1 == j1 {
                // TODO: useless message
                return Err(StructureError::InvalidToken(
                        "self-pairing".to_string(), 
                        "pair list".to_string(), 
                        i1));
            }
            if i1 == 0 || j1 == 0 || i1 > pl.length || j1 > pl.length {
                // TODO: useless message
                return Err(StructureError::InvalidToken(
                        "1-based pair".to_string(), 
                        "pair list".to_string(), 
                        i1));
            }

            let (i, j) = (i1 - 1, j1 - 1);

            if let Some(prev) = table[i] {
                if prev as usize != j {
                    // TODO: useless message
                    return Err(StructureError::UnmatchedOpen(0));
                }
            }
            if let Some(prev) = table[j] {
                if prev as usize != i {
                    // TODO: useless message
                    return Err(StructureError::UnmatchedOpen(0));
                }
            }

            table[i] = Some(j as NAIDX);
            table[j] = Some(i as NAIDX);
        }

        Ok(PairTable(table))
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_list_from_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        let pl = PairList::try_from(&pt).unwrap();

        assert_eq!(pl.length, 6);
        assert_eq!(pl.pairs, vec![Pair(1, 6), Pair(2, 5)]); // 1-based
    }

    #[test]
    fn test_pair_table_from_pair_list_valid() {
        let pl = PairList {
            length: 6,
            pairs: vec![Pair(1, 6), Pair(2, 5)],
        };
        let pt = PairTable::try_from(&pl).unwrap();
        assert_eq!(pt.0, vec![Some(5), Some(4), None, None, Some(1), Some(0)]);
    }

    #[test]
    fn test_pair_table_from_pair_list_out_of_bounds() {
        let pl = PairList {
            length: 4,
            pairs: vec![Pair(1, 5)],
        };
        let err = PairTable::try_from(&pl).unwrap_err();
        assert!(matches!(err, StructureError::InvalidToken(_, _, _)));
    }
}
