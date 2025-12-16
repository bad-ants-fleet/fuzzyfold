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

/// Types of nearest neighbor loops.
///
/// The loops store different combinations of indices (typically u16) that can be used extract
/// relevant Bases from the corresponding sequence. Indices are given in canonical form: (i < j) <
/// (p < q). (There is only one exception to this rule, the ends of JointExterior.)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum NearestNeighborLoop {
    Hairpin {
        closing: (NAIDX, NAIDX), 
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

    /// The "standard" exterior loop, which is the outer-most decomposition of a
    /// single complex 5' to 3' end. For this loop, the index of the 5' end is
    /// smaller than the index of the 3' end.
    Exterior {
        /// The index of the 5' end, followed by the index of the 3' end.
        ends: (NAIDX, NAIDX),
        //NOTE: this list must ALWAYS be in 5'->3' order.
        branches: Vec<(NAIDX, NAIDX)>,
    },

    /// The "joint" exterior loop is one that joins two complexes.
    /// It can be recognized by the fact that the index of the 5' end is 
    /// *larger* than the index of the 3' end.
    /// Evaluation of this loop typically requires rotations out of the 
    /// canonical form, where branches are sorted by index. The first
    /// branch is also the "closing" pair.
    JointExterior { 
        /// The index of the 5' end, followed by the index of the 3' end.
        ends: (NAIDX, NAIDX), //5' end, 3' end.
        //NOTE: this list must ALWAYS be in 5'->3' order.
        branches: Vec<(NAIDX, NAIDX)>,
    },

    /// This is not a NN loop.
    Disconnected {
        left: Box<NearestNeighborLoop>,
        right: Box<NearestNeighborLoop>,
    }
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
            Self::JointExterior { ends: (i, j), branches } => {
                write!(f, "{:<8} [{:>3}, {:>3}], {}", 
                    "Jxterior".red().bold(), i, j,
                    branches.iter()
                    .map(|(i, j)| format!("[{:>3}, {:>3}]", i, j))
                    .collect::<Vec<_>>()
                    .join(", "))
            }
            Self::Disconnected{ .. } => {
                //TODO
                write!(f, "{:<8} ", "Disconnected!".red().bold())
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
            NearestNeighborLoop::JointExterior { branches, .. } => branches.clone(),
            NearestNeighborLoop::Disconnected { .. } => todo!("For now, get this explicitly for the disconnected components!"),
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
            NearestNeighborLoop::JointExterior { ends, branches } => {
                let mut keys = Vec::with_capacity(1 + branches.len());
                keys.push(ends.0 as P1KEY); // NOTE: key only uses 5' end!
                keys.extend(branches.iter().map(|&(i, j)| pack(i, j)));
                keys
            }
            NearestNeighborLoop::Disconnected { .. } => unreachable!("No support for disconnected loop keys."),
        }
    }

    pub fn classify(
        ends: Option<(NAIDX, NAIDX)>, 
        closing: Option<(NAIDX, NAIDX)>, 
        branches: Vec<(NAIDX, NAIDX)>, 
    ) -> Self {
        match closing {
            None => match ends {
                Some((i, j)) if i < j => Self::Exterior { branches, ends: (i, j) },
                Some((i, j)) if j < i => Self::JointExterior { branches, ends: (i, j) },
                _ => panic!("Expected end annotation in exterior loop."),
            },
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
            Self::JointExterior{ branches, .. } => 
                Some(*branches.first().expect("JointExterior must have branches")),
            Self::Exterior { .. } => None,
            Self::Disconnected { .. } => None,
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
            Self::Exterior { ends: (p5, p3), branches } => {
                let mut result = Vec::new();
                let mut start = *p5 as usize;
                for &(p, q) in branches {
                    result.push(start..(p as usize));
                    start = q as usize + 1;
                }
                result.push(start..(*p3 as usize + add1));
                result
            }
            Self::JointExterior { ends: (p5, p3), branches } => {
                debug_assert!(!branches.is_empty());
                // Preprocessing of branches.
                // TODO: more efficient?
                let mut branches = branches.clone();
                branches.rotate_left(1);
                let last = branches.len() - 1;
                let (i, j) = branches[last];
                branches[last] = (j, i);
                while let Some(&(i, _)) = branches.first() {
                    if i > *p3 { break; }
                    branches.rotate_left(1);
                }
                let mut result = Vec::new();
                let mut start = *p5 as usize;
                for (p, q) in branches {
                    result.push(start..(p as usize));
                    start = q as usize + 1;
                }
                result.push(start..(*p3 as usize + add1));
                result
            }
            Self::Disconnected { .. } => todo!("For now, get this explicitly for the disconnected components!"),
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
        debug_assert!(i < j, "Split pair (i,j) must satisfy i < j");
        match self {
            Self::Hairpin { closing: (a, b) } => {
                debug_assert!(*a < i && j < *b, "Pair (i,j) must be within hairpin loop");
                (Self::Interior { closing: (*a, *b), inner: (i, j), },
                 Self::Hairpin { closing: (i, j), })
            }
            
            Self::Interior { closing: (a, b), inner: (p, q) } => {
                debug_assert!(*a < i && j < *b, "Pair (i,j) outside of loop");
                debug_assert!(!(*p < i && j < *q), "Pair (i,j) outside of loop");

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
                (Self::classify(None, Some((*a, *b)), outer_branches),
                 Self::classify(None, Some((i, j)), inner_branches))
            }

            Self::Exterior { ends, branches } => {
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
                (Self::classify(Some(*ends), None, outer_branches),
                 Self::classify(None, Some((i, j)), inner_branches))
            }

            Self::JointExterior { ends: (p5, p3), branches } => {
                let closing = branches.first().expect("JointExterior must have closing pair.");
                let mut outer_branches = vec![(i, j)];
                let mut inner_branches = vec![];
                for &(p, q) in branches.iter().skip(1) {
                    debug_assert!(p < q);
                    if j < p || q < i {
                        outer_branches.push((p, q));
                    } else { 
                        debug_assert!(i < p && q < j);
                        inner_branches.push((p, q));
                    } 
                }
                outer_branches.sort_unstable();
                if i < *p5 && *p3 < j {
                    inner_branches.insert(0, (i, j));
                    (Self::classify(None, Some(*closing), outer_branches),
                    Self::classify(Some((*p5, *p3)), None, inner_branches))
                } else {
                    debug_assert!(j < *p5 || *p3 < i);
                    outer_branches.insert(0, *closing);
                    (Self::classify(Some((*p5, *p3)), None, outer_branches),
                    Self::classify(None, Some((i, j)), inner_branches))
                }
            }

            // TODO: the idea is to only support bimolecular moves.
            Self::Disconnected { .. } => todo!("not sure if we want to join complexes like this."),
        }
    }

    /// Join two loops by reference. 
    /// NOTE: the outer loop must be joined with the inner loop!
    pub fn join_loop(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Disconnected { .. }, _) => { 
                panic!("Cannot join disconnected loops!");
            }

            (_, Self::Disconnected { .. }) => { 
                panic!("Cannot join disconnected loops!");
            }

            (Self::Hairpin { .. }, _) => { 
                panic!("A hairpin cannot be the outer loop!");
            }

            (_, Self::Exterior { .. }) => { 
                panic!("The exterior loop cannot be the inner loop!");
            }

            (Self::JointExterior{ .. }, _) => { 
                panic!("The joint-exterior loop cannot be the outer loop!");
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
 
            (Self::Interior { closing: outer_closing, inner },
             Self::JointExterior { ends, branches }) => {
                assert_eq!(inner, branches.first().unwrap(), 
                    "Cannot join Interior & JointExterior loops!");
                let mut new_branches = Vec::with_capacity(branches.len());
                new_branches.push(*outer_closing);
                new_branches.extend_from_slice(&branches[1..]);
                Self::JointExterior { ends: *ends, branches: new_branches }
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
                let mut new_branches = Vec::with_capacity(
                    outer_branches.len() + inner_branches.len() - 1
                );
                for &x in outer_branches {
                    if x == *inner_closing {
                        new_branches.extend_from_slice(inner_branches);
                    } else {
                        new_branches.push(x);
                    }
                }
                Self::Multibranch { closing: *outer_closing, branches: new_branches }
            },

            (Self::Multibranch { closing: outer_closing, branches: outer_branches },
             Self::JointExterior { ends, branches: inner_branches }
            ) => {
                let inner_closing = inner_branches.first().unwrap();
                let mut new_branches = Vec::with_capacity(
                    outer_branches.len() + inner_branches.len()
                );
                new_branches.push(*outer_closing);
                for &x in outer_branches {
                    if x == *inner_closing {
                        new_branches.extend_from_slice(&inner_branches[1..]);
                    } else {
                        new_branches.push(x);
                    }
                }
                Self::JointExterior { ends: *ends, branches: new_branches }
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
                let mut new_branches =
                    Vec::with_capacity(outer_branches.len() + inner_branches.len() - 1);

                for &x in outer_branches {
                    if x == *inner_closing {
                        new_branches.extend_from_slice(inner_branches);
                    } else {
                        new_branches.push(x);
                    }
                }
                Self::Exterior { ends: *ends, branches: new_branches }
            },

            (Self::Exterior { ends: (lp5, rp3), branches: outer_branches },
             Self::JointExterior { ends: (rp5, lp3), branches: inner_branches }
            ) => {
                let inner_closing = inner_branches.first().unwrap();
                let mut all_branches = Vec::with_capacity(
                    outer_branches.len() + inner_branches.len().saturating_sub(1)
                );

                for &x in outer_branches {
                    if x == *inner_closing {
                        all_branches.extend_from_slice(&inner_branches[1..]);
                    } else {
                        all_branches.push(x);
                    }
                }
                let (left, right): (Vec<_>, Vec<_>) = all_branches
                    .into_iter()
                    .partition(|&(_i, j)| j < *rp5);
 
                let l = Self::Exterior { ends: (*lp5, *lp3), branches: left };
                let r = Self::Exterior { ends: (*rp5, *rp3), branches: right };
                Self::Disconnected { left: Box::new(l), right: Box::new(r)}
            },

        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    macro_rules! assert_closing {
        ($loop:expr, $pair:expr) => {
            assert_eq!($loop.closing(), Some($pair));
        };
    }

    macro_rules! assert_kind {
        ($loop:expr, $pat:pat) => {
            assert!(matches!($loop, $pat), "unexpected loop: {:?}", $loop);
        };
    }


    #[test]
    fn split_hairpin() {
        let lp = NearestNeighborLoop::Hairpin { closing: (2, 10) };

        let (outer, inner) = lp.split_loop(4, 8);

        assert_kind!(outer, NearestNeighborLoop::Interior { .. });
        assert_kind!(inner, NearestNeighborLoop::Hairpin { .. });

        assert_closing!(outer, (2, 10));
        assert_closing!(inner, (4, 8));
    }

    #[test]
    fn split_interior_wrapping() {
        let lp = NearestNeighborLoop::Interior {
            closing: (1, 20),
            inner: (8, 12),
        };

        let (outer, inner) = lp.split_loop(5, 15);

        assert_kind!(outer, NearestNeighborLoop::Interior { .. });
        assert_kind!(inner, NearestNeighborLoop::Interior { .. });

        assert_closing!(outer, (1, 20));
        assert_closing!(inner, (5, 15));
    }

    #[test]
    fn split_interior_left_multibranch() {
        let lp = NearestNeighborLoop::Interior {
            closing: (1, 20),
            inner: (10, 15),
        };

        let (outer, inner) = lp.split_loop(4, 6);

        assert_kind!(outer, NearestNeighborLoop::Multibranch { .. });
        assert_kind!(inner, NearestNeighborLoop::Hairpin { .. });

        assert_closing!(inner, (4, 6));
    }

    #[test]
    fn split_multibranch() {
        let lp = NearestNeighborLoop::Multibranch {
            closing: (1, 30),
            branches: vec![(5, 8), (12, 15), (20, 25)],
        };

        let (outer, inner) = lp.split_loop(10, 18);

        assert_kind!(outer, NearestNeighborLoop::Multibranch { .. });
        assert_kind!(inner, NearestNeighborLoop::Interior { .. });

        assert_closing!(outer, (1, 30));
        assert_closing!(inner, (10, 18));
    }

    #[test]
    fn split_exterior() {
        let lp = NearestNeighborLoop::Exterior {
            ends: (0, 30),
            branches: vec![(5, 8), (12, 15), (20, 25)],
        };

        let (outer, inner) = lp.split_loop(9, 10);

        assert_kind!(outer, NearestNeighborLoop::Exterior { .. });
        assert_kind!(inner, NearestNeighborLoop::Hairpin { .. });

        assert_closing!(inner, (9, 10));
    }

    #[test]
    fn split_joint_exterior_outer() {
        let lp = NearestNeighborLoop::JointExterior {
            ends: (20, 18),
            branches: vec![(1, 28), (5, 8), (12, 15)],
        };

        let (outer, inner) = lp.split_loop(4, 10);

        assert_kind!(outer, NearestNeighborLoop::JointExterior { .. });
        assert_kind!(inner, NearestNeighborLoop::Interior { .. });

        assert_closing!(outer, (1, 28));
        assert_closing!(inner, (4, 10));

        let (outer, inner) = lp.split_loop(4, 25);

        assert_kind!(outer, NearestNeighborLoop::Interior { .. });
        assert_kind!(inner, NearestNeighborLoop::JointExterior { .. });

        assert_closing!(outer, (1, 28));
        assert_closing!(inner, (4, 25));

    }

    #[test]
    fn join_interior_hairpin() {
        let inner = NearestNeighborLoop::Hairpin {
            closing: (3, 7),
        };
        let outer = NearestNeighborLoop::Interior {
            inner: (3, 7),
            closing: (0, 10),
        };
        let joined = outer.join_loop(&inner);
        assert_eq!(
            joined,
            NearestNeighborLoop::Hairpin { closing: (0, 10) }
        );
    }

    #[test]
    fn join_interior_interior() {
        let inner = NearestNeighborLoop::Interior {
            closing: (5, 15),
            inner: (8, 12),
        };
        let outer = NearestNeighborLoop::Interior {
            closing: (0, 20),
            inner: (5, 15),
        };
        let joined = outer.join_loop(&inner);

        assert_eq!(
            joined,
            NearestNeighborLoop::Interior {
                closing: (0, 20),
                inner: (8, 12),
            }
        );
    }

    #[test]
    fn join_interior_multibranch() {
        let outer = NearestNeighborLoop::Interior {
            closing: (0, 20),
            inner: (5, 15),
        };
        let inner = NearestNeighborLoop::Multibranch {
            closing: (5, 15),
            branches: vec![(7, 9), (11, 13)],
        };
        let joined = outer.join_loop(&inner);
        assert_eq!(
            joined,
            NearestNeighborLoop::Multibranch {
                closing: (0, 20),
                branches: vec![(7, 9), (11, 13)],
            }
        );
    }

    #[test]
    fn join_exterior_jointexterior_disconnects() {
        let outer = NearestNeighborLoop::Exterior {
            ends: (0, 25),
            branches: vec![(5, 20)],
        };

        let inner = NearestNeighborLoop::JointExterior {
            ends: (12, 10),
            branches: vec![(5, 20), (6, 10), (15, 19)],
        };

        let joined = outer.join_loop(&inner);

        match joined {
            NearestNeighborLoop::Disconnected { left, right } => {
                assert_eq!(
                    *left,
                    NearestNeighborLoop::Exterior {
                        ends: (0, 10),
                        branches: vec![(6, 10)],
                    }
                );
                assert_eq!(
                    *right,
                    NearestNeighborLoop::Exterior {
                        ends: (12, 25),
                        branches: vec![(15, 19)],
                    }
                );
            }
            other => panic!("expected Disconnected, got {:?}", other),
        }
    }

}

