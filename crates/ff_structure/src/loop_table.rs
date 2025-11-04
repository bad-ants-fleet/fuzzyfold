use std::fmt;
use std::ops::Deref;

use crate::NAIDX;
use crate::PairTable;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopInfo {
    Unpaired { l: NAIDX },
    Paired { o: NAIDX, i: NAIDX }, // outer, inner loop ids
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LoopTable(pub Vec<LoopInfo>);

impl Deref for LoopTable {
    type Target = [LoopInfo];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<&PairTable> for LoopTable {
    fn from(pt: &PairTable) -> Self {
        let n = pt.len();
        let mut table = vec![LoopInfo::Unpaired { l: 0 }; n];
        let mut loop_index: NAIDX = 0;
        let mut mloop: NAIDX = 0;

        let mut stack: Vec<(usize, NAIDX)> = Vec::new(); // (closing_idx, loop_id)
        for i in 0..n {
            match pt[i] {
                None => {
                    table[i] = LoopInfo::Unpaired { l: loop_index };
                }
                Some(j) if (j as usize) > i => {
                    let outer_loop = loop_index;
                    mloop += 1;
                    loop_index = mloop;
                    table[i] = LoopInfo::Paired { o: outer_loop, i: loop_index };
                    stack.push((j as usize, loop_index));
                }
                Some(j) if (j as usize) < i => {
                    let (_, inner_loop) = stack.pop().expect("Expected well formed PairTable, missig opening pair index!");
                    loop_index = stack.last().map(|&(_, l)| l).unwrap_or(0);
                    table[i] = LoopInfo::Paired { o: loop_index, i: inner_loop };
                }
                _ => unreachable!(),
            }
        }
        debug_assert!(stack.is_empty(),
            "Expected well-formed PairTable, missing closing pair index.");
        LoopTable(table)
    }
}

impl fmt::Display for LoopTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut out = Vec::new();
        for info in self.0.iter() {
            let s = match info {
                LoopInfo::Unpaired { l } => format!("{}", l),
                LoopInfo::Paired { o, i } => format!("{}/{}", o, i),
            };
            out.push(s);
        }
        write!(f, "[{}]", out.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_loop_table_valid_structure() {
        // dot-bracket: ((..))
        let pt = PairTable::try_from("((..))").unwrap();
        let lt = LoopTable::from(&pt);

        // Expecting:
        // positions: 0 1 2 3 4 5
        // dotbrack:  ( ( . . ) )
        // loops:    0 1 2   2 1 0
        let expected = [
            LoopInfo::Paired { o: 0, i: 1 },
            LoopInfo::Paired { o: 1, i: 2 },
            LoopInfo::Unpaired { l: 2 },
            LoopInfo::Unpaired { l: 2 },
            LoopInfo::Paired { o: 1, i: 2 },
            LoopInfo::Paired { o: 0, i: 1 },
        ];

        assert_eq!(&lt[..], &expected[..]); // uses Deref
    }

    #[test]
    fn test_loop_table_unpaired_structure() {
        let pt = PairTable::try_from("......").unwrap();
        let lt = LoopTable::from(&pt);

        for info in lt.iter() {
            assert!(matches!(info, LoopInfo::Unpaired { .. }));
        }
    }

    #[test]
    fn test_loop_table_self_pairing_panics() {
        let pt = PairTable::try_from(".").unwrap();
        let lt = LoopTable::from(&pt);
        assert_eq!(lt.len(), 1);
        assert!(matches!(lt[0], LoopInfo::Unpaired { .. }));
    }

    #[test]
    fn test_deref_loop_table_len_indexing() {
        let pt = PairTable::try_from("((..))").unwrap();
        let lt = LoopTable::from(&pt);

        assert_eq!(lt.len(), 6);
        assert!(matches!(lt[2], LoopInfo::Unpaired { .. }));
    }

    #[test]
    fn test_pair_table_to_loop_index_01() {
        use LoopInfo::*;
        let pt = PairTable::try_from(".(((...)).((...))..(.(...)))").unwrap();
        let li = LoopTable( vec![ Unpaired { l: 0 }, 
            Paired { o: 0, i: 1 }, Paired { o: 1, i: 2 }, Paired { o: 2, i: 3 }, 
            Unpaired { l: 3 }, Unpaired { l: 3 }, Unpaired { l: 3 }, 
            Paired { o: 2, i: 3 }, Paired { o: 1, i: 2 },
            Unpaired { l: 1 }, 
            Paired { o: 1, i: 4 }, Paired { o: 4, i: 5 }, 
            Unpaired { l: 5 }, Unpaired { l: 5 }, Unpaired { l: 5 }, 
            Paired { o: 4, i: 5 }, Paired { o: 1, i: 4 }, 
            Unpaired { l: 1 }, Unpaired { l: 1 }, 
            Paired { o: 1, i: 6 }, 
            Unpaired { l: 6 }, 
            Paired { o: 6, i: 7 }, 
            Unpaired { l: 7 }, Unpaired { l: 7 }, Unpaired { l: 7 }, 
            Paired { o: 6, i: 7 }, Paired { o: 1, i: 6 }, Paired { o: 0, i: 1 }
        ]);
        let re = LoopTable::from(&pt);
        assert_eq!(re, li); 
    }

    #[test]
    fn test_pair_table_to_loop_index_02() {
        use LoopInfo::*;
        let pt = PairTable::try_from(".(((...)(...).((.(...))).)).").unwrap();
        let li = LoopTable( vec![ 
                       Unpaired { l: 0 }, 
                       Paired { o: 0, i: 1 },
                       Paired { o: 1, i: 2 },
                       Paired { o: 2, i: 3 },
                       Unpaired { l: 3 },
                       Unpaired { l: 3 }, 
                       Unpaired { l: 3 }, 
                       Paired { o: 2, i: 3 }, 
                       Paired { o: 2, i: 4 }, 
                       Unpaired { l: 4 }, 
                       Unpaired { l: 4 }, 
                       Unpaired { l: 4 },
                       Paired { o: 2, i: 4 }, 
                       Unpaired { l: 2 }, 
                       Paired { o: 2, i: 5 }, 
                       Paired { o: 5, i: 6 }, 
                       Unpaired { l: 6 }, 
                       Paired { o: 6, i: 7 }, 
                       Unpaired { l: 7 }, 
                       Unpaired { l: 7 }, 
                       Unpaired { l: 7 }, 
                       Paired { o: 6, i: 7 }, 
                       Paired { o: 5, i: 6 }, 
                       Paired { o: 2, i: 5 }, 
                       Unpaired { l: 2 }, 
                       Paired { o: 1, i: 2 }, 
                       Paired { o: 0, i: 1 },
                       Unpaired { l: 0 }
        ]);

        let re = LoopTable::from(&pt);
        assert_eq!(re, li); 
    }

    #[test]
    fn test_pair_table_to_loop_index_03() {
        use LoopInfo::*;
        let pt = PairTable::try_from(".(((...)(...))).((((.(...))).)).").unwrap();
        let li = LoopTable( vec![ Unpaired { l: 0 },
                       Paired { o: 0, i: 1 }, Paired { o: 1, i: 2 }, Paired { o: 2, i: 3 }, 
                       Unpaired { l: 3 }, Unpaired { l: 3 }, Unpaired { l: 3 }, 
                       Paired { o: 2, i: 3 }, Paired { o: 2, i: 4 }, 
                       Unpaired { l: 4 }, Unpaired { l: 4 }, Unpaired { l: 4 }, 
                       Paired { o: 2, i: 4 }, Paired { o: 1, i: 2 }, Paired { o: 0, i: 1 }, 
                       Unpaired { l: 0 }, 
                       Paired { o: 0, i: 5 }, Paired { o: 5, i: 6 }, Paired { o: 6, i: 7 }, Paired { o: 7, i: 8 }, 
                       Unpaired { l: 8 }, 
                       Paired { o: 8, i: 9 },
                       Unpaired { l: 9 }, Unpaired { l: 9 }, Unpaired { l: 9 }, 
                       Paired { o: 8, i: 9 }, Paired { o: 7, i: 8 }, Paired { o: 6, i: 7 }, 
                       Unpaired { l: 6 }, Paired { o: 5, i: 6 }, Paired { o: 0, i: 5 }, 
                       Unpaired { l: 0 }
        ]);
        let re = LoopTable::from(&pt);
        assert_eq!(re, li); 
    }

    #[test]
    fn test_loop_table_display() {
        use LoopInfo::*;

        let lt = LoopTable(vec![
            Unpaired { l: 0 },
            Paired { o: 0, i: 1 },
            Paired { o: 1, i: 2 },
            Unpaired { l: 2 },
            Paired { o: 1, i: 2 },
            Paired { o: 0, i: 1 },
        ]);

        let formatted = format!("{}", lt);
        assert_eq!(formatted, "[0, 0/1, 1/2, 2, 1/2, 0/1]");
    }
}

