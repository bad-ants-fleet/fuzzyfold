
use ff_structure::NAIDX;
use ff_structure::PairTable;
use ff_structure::MultiStruct;
use ff_structure::MultiPairTable;
use crate::NearestNeighborLoop;

pub trait LoopDecomposition {
    fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, f: F);

    fn loops(&self) -> Vec<NearestNeighborLoop> {
        let mut out = Vec::new();
        self.for_each_loop(|l| out.push(l.clone()));
        out
    }
}

impl LoopDecomposition for PairTable {
    fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, mut f: F) {
        fn recurse<F: FnMut(&NearestNeighborLoop)>(
            pt: &PairTable,
            closing: Option<(NAIDX, NAIDX)>,
            ends: Option<(NAIDX, NAIDX)>,
            f: &mut F,
        ) {
            let mut branches = Vec::new();

            let (mut p, j) = if let Some((i, j)) = closing {
                (i + 1, j) 
            } else { 
                (0, pt.len() as NAIDX)
            };

            while p < j {
                if let Some(q) = pt[p] {
                    debug_assert!(q > p);
                    branches.push((p, q));
                    recurse(pt, Some((p, q)), None, f);
                    p = q + 1;
                } else {
                    p += 1;
                }
            }
            f(&NearestNeighborLoop::classify(ends, closing, branches));
        }
        recurse(self, None, Some((0, (self.len() - 1) as NAIDX)), &mut f);
    }
}

impl LoopDecomposition for MultiPairTable {
    fn for_each_loop<F: FnMut(&NearestNeighborLoop)>(&self, mut f: F) {
        fn recurse<F: FnMut(&NearestNeighborLoop)>(
            mpt: &MultiPairTable,
            ends: Option<(NAIDX, NAIDX)>,
            closing: Option<(NAIDX, NAIDX)>,
            f: &mut F,
        ) {
            let mut branches = Vec::new();
            let mut closing = closing;
            let mut ends = ends;

            let (mut p, j) = if let Some((i, j)) = closing {
                (i as usize + 1, j as usize) 
            } else { 
                (0, mpt.len())
            };

            while p < j {
                match mpt[p] {
                    MultiStruct::Paired(q) => { 
                        {
                            let p = p as NAIDX;
                            debug_assert!(q > p);
                            branches.push((p, q));
                            recurse(mpt, None, Some((p, q)), f);
                        }
                        p = q as usize + 1;
                    },
                    MultiStruct::Unpaired => p += 1,
                    MultiStruct::StrandBreak => {
                        if let Some((p5, p3)) = ends {
                            let p = p as NAIDX;
                            let inner;
                            (branches, inner) = branches
                                .iter()
                                .copied()
                                .partition(|&(_, k)| k <= p5);
                            f(&NearestNeighborLoop::classify(Some((p5, p - 1)), None, inner));
                            ends = Some((p + 1, p3));
                        } else {
                            let p = p as NAIDX;
                            ends = Some((p + 1, p - 1));
                        }
                        p += 1;
                    },
                }
            }
            if ends.is_some() && let Some((k, l)) = closing {
                branches.insert(0, (k, l));
                closing = None;
            }
            f(&NearestNeighborLoop::classify(ends, closing, branches));
        }
        recurse(self, Some((0 as NAIDX, (self.len() - 1) as NAIDX)), None, &mut f);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::StrandPairTable;

    #[test]
    fn test_decompose_pair_table_edge() {
        let dbn = "."; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![], 
            ends: (0 as NAIDX, 0 as NAIDX), 
        };
        let loops = PairTable::try_from(dbn).expect("valid").loops();
        assert_eq!(loops, vec![eloop]);
    }

    #[test]
    fn test_decompose_pair_table_empty() {
        let dbn = "......."; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![], 
            ends: (0 as NAIDX, 6 as NAIDX), 
        };
        let loops = PairTable::try_from(dbn).expect("valid").loops();
        assert_eq!(loops, vec![eloop]);
    }

    #[test]
    fn test_decompose_multipair_table_empty() {
        let dbn = "......."; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![], 
            ends: (0 as NAIDX, 6 as NAIDX), 
        };
        let spt = StrandPairTable::try_from(dbn).expect("valid");
        let mpt = MultiPairTable::from(&spt);
        assert_eq!(mpt.loops(), vec![eloop]);
    }

    #[test]
    fn test_decompose_multipair_table_single() {
        let dbn = "..(.+..)"; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![(2 as NAIDX, 7 as NAIDX)], 
            ends: (0 as NAIDX, 7 as NAIDX), 
        };
        let jloop = NearestNeighborLoop::JointExterior { 
            branches: vec![(2 as NAIDX, 7 as NAIDX)], 
            ends: (5 as NAIDX, 3 as NAIDX), 
        };

        let spt = StrandPairTable::try_from(dbn).expect("valid");
        let mpt = MultiPairTable::from(&spt);
        assert_eq!(mpt.loops(), vec![jloop, eloop]);
    }

    #[test]
    fn test_decompose_multipair_table_four() {
        let dbn = ".(+(((((.+.))+)))).";
        let mpt = MultiPairTable::try_from(dbn).expect("valid");
        let el1 = NearestNeighborLoop::Exterior { 
            ends: (0 as NAIDX, 18 as NAIDX), 
            branches: vec![(1 as NAIDX, 17 as NAIDX)], 
        };
        let el2 = NearestNeighborLoop::JointExterior { 
            ends: (3 as NAIDX, 1 as NAIDX), 
            branches: vec![(1 as NAIDX, 17 as NAIDX), (3 as NAIDX, 16 as NAIDX)], 
        };
        let el3 = NearestNeighborLoop::JointExterior { 
            ends: (10 as NAIDX, 8 as NAIDX), 
            // -> inv
            branches: vec![(7 as NAIDX, 11 as NAIDX)], 
        };
        let el4 = NearestNeighborLoop::JointExterior { 
            ends: (14 as NAIDX, 12 as NAIDX), 
            branches: vec![(5 as NAIDX, 14 as NAIDX), (6 as NAIDX, 12 as NAIDX)], 
        };
        let loops = mpt.loops();
        println!("{:?}", loops);

        assert!(loops.contains(&el1));
        assert!(loops.contains(&el2));
        assert!(loops.contains(&el3));
        assert!(loops.contains(&el4));

    }


    #[test]
    fn test_loop_ranges() {
        let hloop = NearestNeighborLoop::Hairpin { 
            closing: (1, 5) 
        };
        assert_eq!(vec![2usize, 3, 4], hloop.unpaired_indices());
        assert_eq!(vec![2usize, 3, 4, 5], hloop.inclusive_unpaired_indices());

        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 8) 
        };
        let v: Vec<usize> = Vec::new();
        assert_eq!(v, iloop.unpaired_indices());
        assert_eq!(vec![2usize, 9], iloop.inclusive_unpaired_indices());
 
        let iloop = NearestNeighborLoop::Interior { 
            closing: (1, 9),
            inner: (2, 7) 
        };
        assert_eq!(vec![8usize], iloop.unpaired_indices());
        assert_eq!(vec![2usize, 8, 9], iloop.inclusive_unpaired_indices());
 
        let mloop = NearestNeighborLoop::Multibranch { 
            closing: (1, 15),
            branches: vec![(2, 4), (5, 9)],
        };
        assert_eq!(vec![10usize, 11, 12, 13, 14], mloop.unpaired_indices());
        assert_eq!(vec![2usize, 5, 10, 11, 12, 13, 14, 15], mloop.inclusive_unpaired_indices());

        let eloop = NearestNeighborLoop::Exterior { 
            ends: (0, 15), 
            branches: vec![(1, 5), (6, 11)], 
        };
        assert_eq!(vec![0usize, 12, 13, 14, 15], eloop.unpaired_indices());
        assert_eq!(vec![0usize, 1, 6, 12, 13, 14, 15], eloop.inclusive_unpaired_indices());
 
        let jloop = NearestNeighborLoop::JointExterior { 
            ends: (8, 6), 
            branches: vec![(0, 20), (2, 5), (13, 19)], 
        };
        assert_eq!(vec![8usize, 9, 10, 11, 12, 1, 6], jloop.unpaired_indices());
        assert_eq!(vec![8usize, 9, 10, 11, 12, 13, 20, 1, 2, 6], jloop.inclusive_unpaired_indices());
    }

    #[test]
    fn test_decompose_loops_hairpin() {
        let dbn = ".(...).";
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![(1, 5)], 
            ends: (0, 6), 
        };
        let hloop = NearestNeighborLoop::Hairpin { 
            closing: (1, 5) 
        };
        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
        assert!(loops.len() == 2);
        assert!(loops.contains(&eloop));
        assert!(loops.contains(&hloop));
    }

    #[test]
    fn test_decompose_loops_wild() {
        let dbn = ".(.(((+.(+.)..+))+)).().+...";
        let spt = StrandPairTable::try_from(dbn).expect("valid");
        let loops = MultiPairTable::from(&spt).loops();
        println!("{:?}", loops);

        let el1 = NearestNeighborLoop::Exterior { 
            ends: (0, 23), 
            branches: vec![(1, 19), (21, 22)], 
        };
        assert!(loops.contains(&el1));

        let il1 = NearestNeighborLoop::Interior {
            closing: (1, 19), 
            inner: (3, 18) 
        };
        assert!(loops.contains(&il1));

        let el2 = NearestNeighborLoop::Exterior { 
            ends: (25, 27), 
            branches: vec![],
        };
        assert!(loops.contains(&el2));

        let jl1 = NearestNeighborLoop::JointExterior { 
            ends: (10, 8), 
            branches: vec![(8, 11)],
        };
        assert!(loops.contains(&jl1));
            
        let el3 = NearestNeighborLoop::Exterior { 
            ends: (7, 13), 
            branches: vec![(8, 11)],
        };
        assert!(loops.contains(&el3));

        let jl2 = NearestNeighborLoop::JointExterior { 
            ends: (18, 16), 
            branches: vec![(3, 18), (4, 16)],
        };
        assert!(loops.contains(&jl2));

        let jl3 = NearestNeighborLoop::JointExterior { 
            ends: (15, 5), 
            branches: vec![(5, 15)],
        };
        assert!(loops.contains(&jl3));

    }
}

