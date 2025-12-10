//! The nearest neighbor loop representations. 
//!
//! We may update this to use ff_structure::Pair in the future,
//! so that we get (NAIDX, NAIDX) -> P1KEY conversions out of the box.
//!

use std::fmt;
use std::ops::Range;
use colored::*;

use ff_structure::NAIDX;
use ff_structure::P1KEY;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum NearestNeighborLoop {
    Hairpin {
        closing: (NAIDX, NAIDX), // (i, j)
    },
    Interior {
        closing: (NAIDX, NAIDX),
        inner: (NAIDX, NAIDX),
    },
    Multibranch {
        closing: (NAIDX, NAIDX),
        //NOTE: this list must ALWAYS be in 5'->3' order.
        branches: Vec<(NAIDX, NAIDX)>,
    },
    Exterior {
        //NOTE: this list must ALWAYS be in 5'->3' order.
        ends: (NAIDX, NAIDX),
        branches: Vec<(NAIDX, NAIDX)>,
    },
}

impl fmt::Display for NearestNeighborLoop {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Hairpin { closing: (i, j) } => {
                write!(f, "{:<8} ({:>3}, {:>3})", 
                    "Hairpin".cyan(), i, j)
            }
            Self::Interior { closing: (i, j), inner: (p, q) } => {
                write!(f, "{:<8} ({:>3}, {:>3}), ({:>3}, {:>3})", 
                    "Interior".cyan(), i, j, p, q)
            }
            Self::Multibranch { closing: (i, j), branches } => {
                write!(f, "{:<8} ({:>3}, {:>3}), {}", 
                    "Multibr.".cyan().bold(), i, j, 
                    branches.iter()
                    .map(|(i, j)| format!("[{:>3}, {:>3}]", i, j))
                    .collect::<Vec<_>>()
                    .join(", "))
            }
            Self::Exterior { ends: (i, j), branches } => {
                write!(f, "{:<8} [{:>3}, {:>3}], {}", 
                    "Exterior".cyan().bold(), i, j,
                    branches.iter()
                    .map(|(i, j)| format!("[{:>3}, {:>3}]", i, j))
                    .collect::<Vec<_>>()
                    .join(", "))
            }
        }
    }
}

impl NearestNeighborLoop {

    /// Return all base pairs (closing, inner, and/or branches) adjacent to this loop.
    pub fn pairs(&self) -> Vec<(NAIDX, NAIDX)> {
        match self {
            NearestNeighborLoop::Hairpin { closing } => vec![*closing],
            NearestNeighborLoop::Interior { closing, inner } => vec![*closing, *inner],
            NearestNeighborLoop::Multibranch { closing, branches } => {
                let mut pairs = Vec::with_capacity(1 + branches.len());
                pairs.push(*closing);
                pairs.extend(branches.iter().cloned());
                pairs
            }
            NearestNeighborLoop::Exterior { branches, .. } => branches.clone(),
        }
    }

    /// Return all base pairs (closing, inner, and/or branches) as packed P1KEY.
    /// Exterior loops use 0 as a sentinel closing key.
    pub fn loop_key(&self) -> Vec<P1KEY> {
        fn pack(i: NAIDX, j: NAIDX) -> P1KEY {
            ((i as P1KEY) << NAIDX::BITS) | (j as P1KEY)
        }

        match self {
            NearestNeighborLoop::Hairpin { closing } => {
                vec![pack(closing.0, closing.1)]
            }
            NearestNeighborLoop::Interior { closing, inner } => {
                vec![
                    pack(closing.0, closing.1),
                    pack(inner.0, inner.1),
                ]
            }
            NearestNeighborLoop::Multibranch { closing, branches } => {
                let mut keys = Vec::with_capacity(1 + branches.len());
                keys.push(pack(closing.0, closing.1));
                keys.extend(branches.iter().map(|&(i, j)| pack(i, j)));
                keys
            }
            NearestNeighborLoop::Exterior { ends, branches } => {
                let mut keys = Vec::with_capacity(1 + branches.len());
                keys.push(ends.0 as P1KEY); // NOTE: key only uses 5' end!
                keys.extend(branches.iter().map(|&(i, j)| pack(i, j)));
                keys
            }
        }
    }

    pub fn classify(
        ends: Option<(NAIDX, NAIDX)>, 
        closing: Option<(NAIDX, NAIDX)>, 
        branches: Vec<(NAIDX, NAIDX)>, 
    ) -> Self {
        match closing {
            None => Self::Exterior { branches, ends: ends.expect("check") },
            Some((i, j)) => match branches.len() {
                0 => Self::Hairpin { closing: (i, j) },
                1 => Self::Interior { closing: (i, j), inner: branches[0] },
                _ => Self::Multibranch { closing: (i, j), branches },
            },
        }
    }

    pub fn closing(&self) -> Option<(NAIDX, NAIDX)> {
        match self { Self::Hairpin { closing }
            | Self::Interior { closing, .. }
            | Self::Multibranch { closing, .. } => Some(*closing),
            Self::Exterior { .. } => None,
        }
    }

    /// Returns a list of ranges for unpaired nucleotides.
    ///
    /// NOTE: add1 feels a bit like a hack. It is used by the calling function
    /// to ensure that the exterior loop always stops at the last unpaired
    /// nucleotide.
    fn unpaired_ranges(&self, add1:  usize) -> Vec<Range<usize>> {
        match self {
            Self::Hairpin { closing: (i, j) } => vec![
                (*i as usize + 1)..(*j as usize)],
            Self::Interior { closing: (i, j),  inner: (p, q) } => vec![
                (*i as usize + 1)..(*p as usize), 
                (*q as usize + 1)..(*j as usize)
            ],
            Self::Multibranch { closing: (i, j), branches } => {
                let mut result = vec![];
                let mut start = *i as usize;
                for &(p, q) in branches {
                    result.push((start + 1)..(p as usize));
                    start = q as usize;
                }
                result.push((start+1)..(*j as usize));
                result
            }
            Self::Exterior { ends: (i, j), branches } => {
                let mut result = Vec::new();
                let mut start = *i as usize;
                for &(p, q) in branches {
                    result.push(start..(p as usize));
                    start = q as usize + 1;
                }
                result.push(start..(*j as usize + add1));
                result
            }
        }
    }

    pub fn unpaired_indices(&self) -> Vec<usize> {
        self.unpaired_ranges(1)
            .into_iter()
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect()
    }

    /// Returns all sequence indices that should point to this loop.
    pub fn inclusive_unpaired_indices(&self) -> Vec<usize> {
        self.unpaired_ranges(0)
            .into_iter()
            .map(|r| r.start..=r.end)
            .flat_map(|r| r.collect::<Vec<_>>())
            .collect()
    }

    /// Split the given loop into two new loops at the indices i,j
    /// NOTE: Returns (outer, inner)
    pub fn split_loop(&self, i: NAIDX, j: NAIDX) -> (Self, Self) {
        assert!(i < j, "Split pair (i,j) must satisfy i < j");
        match self {
            Self::Hairpin { closing: (a, b) } => {
                assert!(*a < i && j < *b, "Pair (i,j) must be within hairpin loop");
                (Self::Interior { closing: (*a, *b), inner: (i, j), },
                 Self::Hairpin { closing: (i, j), })
            }
            
            Self::Interior { closing: (a, b), inner: (p, q) } => {
                assert!(*a < i && j < *b, "Pair (i,j) outside of loop");
                assert!(!(*p < i && j < *q), "Pair (i,j) outside of loop");

                if i < *p && *q < j {
                    (Self::Interior { closing: (*a, *b), inner: (i, j) },
                     Self::Interior { closing: (i, j), inner: (*p, *q) })
                } else if j < *p {
                    (Self::Multibranch { closing: (*a, *b), branches: vec![(i, j), (*p, *q)] },
                     Self::Hairpin { closing: (i, j) })
                } else if *q < i {
                    (Self::Multibranch { closing: (*a, *b), branches: vec![(*p, *q), (i, j)] },
                     Self::Hairpin { closing: (i, j) })
                } else {
                    panic!("that really should not happen.");
                }
            }

            Self::Multibranch { closing: (a, b), branches } => {
                debug_assert!(*a < i && j < *b, "Pair (i,j) outside loop");

                let mut outer_branches = vec![(i, j)];
                let mut inner_branches = vec![];
    
                for &(p, q) in branches {
                    debug_assert!(p < q);
                    if j < p || q < i {
                        outer_branches.push((p, q));
                    } else { 
                        debug_assert!(i < p && q < j);
                        inner_branches.push((p, q));
                    } 
                }

                outer_branches.sort_unstable();
                inner_branches.sort_unstable();
                (Self::classify(None, Some((*a, *b)), outer_branches),
                 Self::classify(None, Some((i, j)), inner_branches))
            }

            Self::Exterior { ends, branches } => {
                let mut outer_branches = vec![];
                let mut inner_branches = vec![];
    
                outer_branches.push((i, j));
                for &(p, q) in branches {
                    debug_assert!(p < q);
                    if j < p || q < i {
                        outer_branches.push((p, q));
                    } else { 
                        debug_assert!(i < p && q < j);
                        inner_branches.push((p, q));
                    } 
                }

                outer_branches.sort_unstable();
                inner_branches.sort_unstable();
                (Self::classify(Some(*ends), None, outer_branches),
                 Self::classify(None, Some((i, j)), inner_branches))
            }
        }
    }

    /// Join two loops by reference. 
    /// NOTE: the outer loop must be joined with the inner loop!
    pub fn join_loop(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Hairpin { .. }, _) => { 
                panic!("A hairpin cannot be the outer loop!");
            }

            (_, Self::Exterior { .. }) => { 
                panic!("Multi-stranded moves are not supported.");
            }

            (Self::Interior { closing: outer_closing, inner },
             Self::Hairpin { closing: inner_closing }) => {
                assert_eq!(inner, inner_closing, "Cannot join interior & haipin loops!");
                Self::Hairpin { closing: *outer_closing } 
            }

            (Self::Interior { closing: outer_closing, inner: outer_inner },
             Self::Interior { closing: inner_closing, inner: inner_inner })  => {
                assert_eq!(outer_inner, inner_closing, "Cannot join interior & interior loops!");
                Self::Interior { closing: *outer_closing, inner: *inner_inner }
            }
 
            (Self::Interior { closing: outer_closing, inner },
             Self::Multibranch { closing: inner_closing, branches }) => {
                assert_eq!(inner, inner_closing, "Cannot join interior & multibranch loops!");
                Self::Multibranch { closing: *outer_closing, branches: branches.clone() }
            },

            (Self::Multibranch { closing: outer_closing, branches },
             Self::Hairpin { closing: inner_closing }
            ) => {
                let new_branches: Vec<_> = branches.iter().cloned()
                    .filter(|b| b != inner_closing)
                    .collect();
                assert_eq!(branches.len(), new_branches.len() + 1, "Cannot join multibranch & hairpin loops!");
                Self::classify(None, Some(*outer_closing), new_branches)
            }
  
            (Self::Multibranch { closing: outer_closing, branches },
             Self::Interior { closing: inner_closing, inner }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .map(|x| if x == *inner_closing { *inner } else { x })
                    .collect();
                Self::Multibranch { closing: *outer_closing, branches: new_branches }
            },
 
            (Self::Multibranch { closing: outer_closing, branches: outer_branches },
             Self::Multibranch { closing: inner_closing, branches: inner_branches }
            ) => {
                let new_branches: Vec<_>  = outer_branches.iter().cloned()
                    .flat_map(|x| if x == *inner_closing { inner_branches.clone()  } else { vec![x] })
                    .collect();
                Self::Multibranch { closing: *outer_closing, branches: new_branches }
            },
            
            (Self::Exterior { ends, branches },
             Self::Hairpin { closing: inner_closing }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .filter(|b| b != inner_closing)
                    .collect();
                Self::Exterior { ends: *ends, branches: new_branches }
            }

            (Self::Exterior { ends, branches },
             Self::Interior { closing: inner_closing, inner }
            ) => {
                let new_branches: Vec<_>  = branches.iter().cloned()
                    .map(|x| if x == *inner_closing { *inner } else { x })
                    .collect();
                Self::Exterior { ends: *ends, branches: new_branches }
            },
              
            (Self::Exterior { ends, branches: outer_branches },
             Self::Multibranch { closing: inner_closing, branches: inner_branches }
            ) => {
                let new_branches: Vec<_>  = outer_branches.iter().cloned()
                    .flat_map(|x| if x == *inner_closing { inner_branches.clone()  } else { vec![x] })
                    .collect();
                Self::Exterior { ends: *ends, branches: new_branches }
            },
        }
    }
}
