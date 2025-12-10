
use ff_structure::NAIDX;
use ff_structure::PairTable;
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
                (i as usize + 1, j as usize) 
            } else { 
                (0, pt.len())
            };

            while p < j {
                if let Some(q) = pt[p] {
                    debug_assert!(q > p as NAIDX);
                    branches.push((p as NAIDX, q));
                    recurse(pt, Some((p as NAIDX, q)), None, f);
                    p = q as usize + 1;
                } else {
                    p += 1;
                }
            }
            f(&NearestNeighborLoop::classify(ends, closing, branches));
        }
        recurse(self, None, Some((0 as NAIDX, (self.len() - 1) as NAIDX)), &mut f);
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decompose_loops_empty() {
        let dbn = "......."; // all unpaired → exterior loop only
        let eloop = NearestNeighborLoop::Exterior { 
            branches: vec![], 
            ends: (0 as NAIDX, 6 as NAIDX), 
        };

        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
        assert_eq!(loops, vec![eloop]);
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
        let dbn = ".(.((...))()((()))).((...()))";
        let loops = PairTable::try_from(dbn).expect("valid").loops();
        println!("{:?}", loops);
    }
}

