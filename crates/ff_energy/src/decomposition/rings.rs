//! The ring decomposition: loops nested along the backbone. From a paired
//! base, the walk continues after its partner. Without crossing pairs, the
//! loops are found pair by pair, between its two bases (`loops_closed_by`);
//! with crossing pairs, all loops are nested, and pseudoknots are the loops
//! on cycles of the loop graph (see `pseudoknots.rs`).

use std::fmt;

use ff_structure::MultiPairTable;
use ff_structure::MultiStruct;

use crate::Pseudoknot;
use crate::Gap;
use crate::NestedLoop;
use std::mem::take;

use super::nested_loops::Pairs;
use super::pseudoknots::trace_pseudoknots;

/// A part of a ring decomposition: a loop, or a pseudoknot.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Ring {
    Loop(NestedLoop),
    Pseudoknot(Pseudoknot),
}

impl fmt::Display for Ring {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Loop(l) => fmt::Display::fmt(l, f),
            Self::Pseudoknot(pk) => fmt::Display::fmt(pk, f),
        }
    }
}

/// A structure split into loops and pseudoknots: what an
/// `EnergyModel` evaluates.
pub trait RingDecomposition {
    /// Call `f` on every loop and pseudoknot. Without crossing pairs, a loop
    /// after the loops it encloses, the loops of the exterior (open loops
    /// and circles) last; with crossing pairs, the loops on no cycle, then
    /// the pseudoknots.
    fn for_each_ring<F: FnMut(&Ring)>(&self, f: F);

    /// All loops and pseudoknots, in the order of [`Self::for_each_ring`].
    fn rings(&self) -> Vec<Ring> {
        let mut out = Vec::new();
        self.for_each_ring(|r| out.push(r.clone()));
        out
    }
}

impl RingDecomposition for MultiPairTable {
    fn for_each_ring<F: FnMut(&Ring)>(&self, mut f: F) {
        if !self.is_nested() {
            return trace_pseudoknots(self, &mut f);
        }
        trace_nested_loops(self, |l| f(&Ring::Loop(l)));
    }
}

/// The loops of a structure without crossing pairs (see
/// `MultiPairTable::is_nested`): those closed by each pair, 5' to 3' by its
/// 3' end, then those of the exterior.
fn trace_nested_loops<F: FnMut(NestedLoop)>(mpt: &MultiPairTable, mut f: F) {
    debug_assert!(mpt.is_nested(), "Crossing pairs: see trace_pseudoknots.");
    for k in 0..mpt.len() {
        if let MultiStruct::Paired(p) = mpt.get(k) && p < k {
            loops_closed_by(mpt, Some((p, k)), &mut f);
        }
    }
    loops_closed_by(mpt, None, &mut f);
}

/// The loops between i and j, outside the pairs they enclose: the loop
/// closed by (i, j), and with strand ends, the loops they split off. A
/// strand break ends an open loop, a strand wrap closes a circle. With
/// `None`, the loops of the exterior.
fn loops_closed_by<F: FnMut(NestedLoop)>(mpt: &MultiPairTable, closing: Option<(usize, usize)>, f: &mut F) {
    let (lo, hi) = closing.map_or((0, mpt.len()), |(i, j)| (i + 1, j));
    // The pairs since the last strand end, at `five`.
    let (mut five, mut pairs) = (lo, Pairs::None);
    // The pairs and the 3' end before the first strand break, joined to the
    // last strand by the closing pair.
    let mut head: Option<(Pairs, usize)> = None;
    let mut k = lo;
    while k < hi {
        match mpt.get(k) {
            MultiStruct::Paired(q) => {
                pairs.push((k, q));
                k = q;
            }
            MultiStruct::StrandWrap => {
                f(NestedLoop::new(take(&mut pairs), Gap::Wrap { last: k - 1, first: five }));
                five = k + 1;
            }
            MultiStruct::StrandBreak if closing.is_some() && head.is_none() => {
                head = Some((take(&mut pairs), k - 1));
                five = k + 1;
            }
            MultiStruct::StrandBreak => {
                f(NestedLoop::new(take(&mut pairs), Gap::Open { five, three: k - 1 }));
                five = k + 1;
            }
            MultiStruct::Unpaired => (),
        }
        k += 1;
    }
    match (closing, head) {
        (Some((i, j)), None) => {
            pairs.push((j, i));
            f(NestedLoop::new(pairs, Gap::None));
        }
        (Some((i, j)), Some((head, three))) => {
            pairs.push((j, i));
            pairs.extend(&head);
            f(NestedLoop::new(pairs, Gap::Open { five, three }));
        }
        // Without a linear strand at the 3' end, nothing is left.
        (None, _) if five < hi => f(NestedLoop::new(pairs, Gap::Open { five, three: hi - 1 })),
        (None, _) => (),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lp(pairs: &[(usize, usize)], gap: Gap) -> Ring {
        Ring::Loop(NestedLoop::new(Pairs::from(pairs.to_vec()), gap))
    }

    #[test]
    fn loops_closed_by_a_pair() {
        let closed_by = |db: &str, closing: Option<(usize, usize)>| {
            let mut out = Vec::new();
            loops_closed_by(&MultiPairTable::try_from(db).unwrap(), closing, &mut |l| out.push(Ring::Loop(l)));
            out
        };
        assert_eq!(closed_by("(...)", Some((0, 4))), [lp(&[(4, 0)], Gap::None)]);
        // A strand break: the closing pair joins the last strand to the first.
        assert_eq!(closed_by("((..+..))", Some((1, 7))), [lp(&[(7, 1)], Gap::Open { five: 5, three: 3 })]);
        // Two strand breaks split off the strand between them.
        assert_eq!(closed_by("(.+.(.).+.)", Some((0, 10))), [
            lp(&[(4, 6)], Gap::Open { five: 3, three: 7 }),
            lp(&[(10, 0)], Gap::Open { five: 9, three: 1 }),
        ]);
        // A circle inside: closed by its wrap, one strand end around it.
        assert_eq!(closed_by("((..+....^..))", Some((1, 12))), [
            lp(&[], Gap::Wrap { last: 8, first: 5 }),
            lp(&[(12, 1)], Gap::Open { five: 10, three: 3 }),
        ]);
        // The exterior: a circle before the linear strand is not part of it.
        assert_eq!(closed_by("..^.(.)", None), [
            lp(&[], Gap::Wrap { last: 1, first: 0 }),
            lp(&[(4, 6)], Gap::Open { five: 3, three: 6 }),
        ]);
    }

    #[test]
    fn loops_of_nested_structures() {
        for (db, rings) in [
            (".", vec![
                lp(&[], Gap::Open { five: 0, three: 0 }),
            ]),
            (".(...).", vec![
                lp(&[(5, 1)], Gap::None),
                lp(&[(1, 5)], Gap::Open { five: 0, three: 6 }),
            ]),
            // A strand break inside a pair: the loop through it is open.
            ("..(.+..)", vec![
                lp(&[(7, 2)], Gap::Open { five: 5, three: 3 }),
                lp(&[(2, 7)], Gap::Open { five: 0, three: 7 }),
            ]),
            (".(+(((((.+.))+)))).", vec![
                lp(&[(11, 7)], Gap::Open { five: 10, three: 8 }),
                lp(&[(7, 11), (12, 6)], Gap::None),
                lp(&[(14, 5), (6, 12)], Gap::Open { five: 14, three: 12 }),
                lp(&[(5, 14), (15, 4)], Gap::None),
                lp(&[(4, 15), (16, 3)], Gap::None),
                lp(&[(3, 16), (17, 1)], Gap::Open { five: 3, three: 1 }),
                lp(&[(1, 17)], Gap::Open { five: 0, three: 18 }),
            ]),
            // Two strand breaks inside a pair split off an open loop.
            (".(.(((+.(+.)..+))+)).().+...", vec![
                lp(&[(11, 8)], Gap::Open { five: 10, three: 8 }),
                lp(&[(8, 11)], Gap::Open { five: 7, three: 13 }),
                lp(&[(15, 5)], Gap::Open { five: 15, three: 5 }),
                lp(&[(5, 15), (16, 4)], Gap::None),
                lp(&[(18, 3), (4, 16)], Gap::Open { five: 18, three: 16 }),
                lp(&[(3, 18), (19, 1)], Gap::None),
                lp(&[(22, 21)], Gap::None),
                lp(&[(1, 19), (21, 22)], Gap::Open { five: 0, three: 23 }),
                lp(&[], Gap::Open { five: 25, three: 27 }),
            ]),
            // The loop through the wrap of a circle, also without pairs.
            ("....^", vec![
                lp(&[], Gap::Wrap { last: 3, first: 0 }),
            ]),
            (".(..)(..).^", vec![
                lp(&[(4, 1)], Gap::None),
                lp(&[(8, 5)], Gap::None),
                lp(&[(1, 4), (5, 8)], Gap::Wrap { last: 9, first: 0 }),
            ]),
            // A circle between two linear strands: closed by its wrap, and one
            // strand end in the loop around it.
            ("((..+....^..))", vec![
                lp(&[], Gap::Wrap { last: 8, first: 5 }),
                lp(&[(12, 1)], Gap::Open { five: 10, three: 3 }),
                lp(&[(1, 12), (13, 0)], Gap::None),
                lp(&[(0, 13)], Gap::Open { five: 0, three: 13 }),
            ]),
            // Circles before and after the linear strand: not in its exterior loop.
            ("..^.(...).+..(...)^", vec![
                lp(&[(8, 4)], Gap::None),
                lp(&[(17, 13)], Gap::None),
                lp(&[], Gap::Wrap { last: 1, first: 0 }),
                lp(&[(4, 8)], Gap::Open { five: 3, three: 9 }),
                lp(&[(13, 17)], Gap::Wrap { last: 17, first: 11 }),
            ]),
        ] {
            assert_eq!(MultiPairTable::try_from(db).unwrap().rings(), rings, "{db}");
        }
    }

    /// The loops of the walk between paired bases and of the general tracer
    /// (see `pseudoknots.rs`) are the same, up to their order.
    fn check_walk(db: &str) {
        let mpt = MultiPairTable::try_from(db).unwrap();
        assert!(mpt.is_nested(), "{db}");
        let mut walked = Vec::new();
        trace_nested_loops(&mpt, |l| walked.push(format!("{l:?}")));
        let mut nested = Vec::new();
        trace_pseudoknots(&mpt, &mut |r: &Ring| match r {
            Ring::Loop(l) => nested.push(format!("{l:?}")),
            Ring::Pseudoknot(pk) => panic!("{db}: {pk}"),
        });
        walked.sort();
        nested.sort();
        assert_eq!(walked, nested, "{db}");
    }

    #[test]
    fn walk_equals_general_tracer() {
        for db in [
            ".", ".(...).", "..(.+..)", ".(+(((((.+.))+)))).", ".(.(((+.(+.)..+))+)).().+...",
            "((((...))))......^", "....^", "((..+....^..))", "(..+.(...).^.(.).^..)", "..^...^",
            ".(...).^((...))", "(...).+..^", "..^.(...).+..(...)^", "((..+..((..+..))..+..))",
        ] {
            check_walk(db);
        }
        let quarters = |s: &str| {
            let n = s.len();
            format!("{}&{}&{}&{}", &s[..n / 4], &s[n / 4..n / 2], &s[n / 2..3 * n / 4], &s[3 * n / 4..])
        };
        for len in [50, 100, 250] {
            let path = format!("{}/benches/data/benchmark_random_structures_len{len}.vrna", env!("CARGO_MANIFEST_DIR"));
            for db in std::fs::read_to_string(path).unwrap().lines().skip(2).step_by(3) {
                for db in [db.to_string(), quarters(db), format!("{db}^"), format!("{}^{}", &db[..len / 2], &db[len / 2..])] {
                    if MultiPairTable::try_from(db.as_str()).unwrap().is_nested() {
                        check_walk(&db);
                    }
                }
            }
        }
    }

}
