//! Loops nested along the backbone (see `RingDecomposition`).

use std::borrow::Cow;
use std::fmt;
use std::ops::Range;

use crate::Base;

type Pair = (usize, usize);

/// The pairs of a loop: none, one or two inline, more in a `Vec`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub(crate) enum Pairs {
    #[default]
    None,
    One((usize, usize)),
    Two([(usize, usize); 2]),
    Many(Vec<(usize, usize)>),
}

impl From<Vec<(usize, usize)>> for Pairs {
    fn from(pairs: Vec<(usize, usize)>) -> Self {
        match pairs.as_slice() {
            [] => Pairs::None,
            [p] => Pairs::One(*p),
            [p, q] => Pairs::Two([*p, *q]),
            _ => Pairs::Many(pairs),
        }
    }
}

impl Pairs {
    pub(crate) fn as_slice(&self) -> &[(usize, usize)] {
        match self {
            Pairs::None => &[],
            Pairs::One(p) => std::slice::from_ref(p),
            Pairs::Two(p) => p,
            Pairs::Many(p) => p,
        }
    }

    fn as_mut_slice(&mut self) -> &mut [(usize, usize)] {
        match self {
            Pairs::None => &mut [],
            Pairs::One(p) => std::slice::from_mut(p),
            Pairs::Two(p) => p,
            Pairs::Many(p) => p,
        }
    }

    /// Rotate the pairs `k` to the left.
    pub(crate) fn rotate_left(&mut self, k: usize) {
        self.as_mut_slice().rotate_left(k);
    }

    /// The pairs of all parts, in order.
    fn concat(parts: &[&[(usize, usize)]]) -> Self {
        let n = parts.iter().map(|p| p.len()).sum();
        let mut all = parts.iter().flat_map(|p| p.iter().copied());
        match n {
            0 => Pairs::None,
            1 => Pairs::One(all.next().unwrap()),
            2 => Pairs::Two([all.next().unwrap(), all.next().unwrap()]),
            _ => {
                let mut v = Vec::with_capacity(n);
                parts.iter().for_each(|p| v.extend_from_slice(p));
                Pairs::Many(v)
            }
        }
    }

    pub(crate) fn push(&mut self, pair: (usize, usize)) {
        match self {
            Pairs::Many(v) => v.push(pair),
            Pairs::None => *self = Pairs::One(pair),
            Pairs::One(p) => *self = Pairs::Two([*p, pair]),
            Pairs::Two([p, q]) => *self = Pairs::Many(vec![*p, *q, pair]),
        }
    }

    /// Append the pairs of `other`.
    pub(crate) fn extend(&mut self, other: &Pairs) {
        let all = Pairs::concat(&[self.as_slice(), other.as_slice()]);
        *self = all;
    }
}

/// A stretch of a loop between two helix ends (see `NestedLoop::stretches`):
/// its ends, the paired bases (or the strand ends of an open loop or a
/// circle), and its unpaired bases. The bases that point to the loop, for
/// the lookup of the move sets, are `unpaired.start..=ends.1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Stretch {
    pub ends: (usize, usize),
    pub unpaired: Range<usize>,
}

/// Where the walk around a loop leaves the backbone, if anywhere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Gap {
    /// A closed loop.
    None,
    /// An open loop, from the 5' end of a strand to the 3' end of a strand.
    Open { five: usize, three: usize },
    /// The loop through the strand wrap of a circular strand: its last and
    /// its first base.
    Wrap { last: usize, first: usize },
}

/// A loop nested along the backbone: its pairs in loop order, each (a, b)
/// as seen from the loop, where the walk arrives at a and continues at b.
/// A closed loop lists its closing pair (j, i) last; an open loop starts at
/// its 5' end, a loop through the wrap at the first base of the strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NestedLoop {
    pairs: Pairs,
    gap: Gap,
}

impl NestedLoop {
    pub(crate) fn new(pairs: Pairs, gap: Gap) -> Self {
        Self { pairs, gap }
    }

    /// The pairs in loop order, each (a, b) as seen from the loop.
    pub fn pairs(&self) -> &[(usize, usize)] {
        self.pairs.as_slice()
    }

    /// Where the walk around the loop leaves the backbone.
    pub fn gap(&self) -> Gap {
        self.gap
    }

    /// The segments of the loop in loop order, each from the base where the
    /// walk continues to the base where it arrives: for a closed loop and
    /// the loop through a wrap, from the last pair to the first, then
    /// between consecutive pairs; for an open loop, from its 5' end to the
    /// first pair, between consecutive pairs, and from the last pair to its
    /// 3' end. A segment that passes a wrap is a copy; from a pair back to
    /// the same base, it is the whole circle. Without pairs, the one
    /// segment of the loop.
    pub fn segments<'a>(&'a self, sequence: &'a [Base]) -> impl Iterator<Item = Cow<'a, [Base]>> + 'a {
        // From x to y, through the wrap of a circle unless `forward`.
        let segment = move |(x, y, forward): (usize, usize, bool)| -> Cow<'a, [Base]> {
            if forward {
                return Cow::Borrowed(&sequence[x..=y]);
            }
            let end = sequence[x..].iter().position(|&b| b == Base::SW).map_or(sequence.len(), |e| x + e);
            let first = sequence[..y].iter().rposition(|&b| b == Base::SB || b == Base::SW).map_or(0, |e| e + 1);
            Cow::Owned(sequence[x..end].iter().chain(&sequence[first..=y]).copied().collect())
        };
        let pairs = self.pairs();
        // The segment before the first pair, and for an open loop the one after the last.
        let (head, tail) = match (self.gap, pairs) {
            (Gap::Open { five, three }, []) => ((five, three, true), None),
            (Gap::Wrap { last, first }, []) => ((first, last, true), None),
            (Gap::Open { five, three }, [p, .., q] | [p @ q]) =>
                ((five, p.0, five <= p.0), Some((q.1, three, q.1 <= three))),
            (Gap::None | Gap::Wrap { .. }, [p, .., q] | [p @ q]) => ((q.1, p.0, q.1 < p.0), None),
            (Gap::None, []) => unreachable!("A closed loop has a closing pair."),
        };
        let between = pairs.windows(2).map(|w| (w[0].1, w[1].0, w[0].1 < w[1].0));
        std::iter::once(head).chain(between).chain(tail).map(segment)
    }
}

/// By loop type (hairpin, interior, multiloop, exterior, joint exterior,
/// circular), the pairs (i, j) with i < j; a joint exterior loop (from a 5'
/// end after its 3' end) lists them in loop order. A loop that meets its
/// pairs out of order (of a pseudoknot) lists them in loop order as seen
/// from the loop.
impl fmt::Display for NestedLoop {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let canonical: Vec<(usize, usize)> = self.pairs().iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        // Branches in index order, and only the closing pair of a closed loop seen from inside.
        let body = match self.gap { Gap::None => &self.pairs()[..self.pairs().len().saturating_sub(1)], _ => self.pairs() };
        let ordered = body.windows(2).all(|w| w[0].1 < w[1].0) && body.iter().all(|&(a, b)| a < b);
        let joint = matches!(self.gap, Gap::Open { five, three } if five > three);
        if !ordered && !joint {
            let pairs: Vec<String> = self.pairs().iter().map(|(a, b)| format!("({a:>3}, {b:>3})")).collect();
            return write!(f, "{:<8} {}", "Pseudol.", pairs.join(", "));
        }
        match (self.gap, self.pairs()) {
            (Gap::None, [(j, i)]) => write!(f, "{:<8} ({:>3}, {:>3})", "Hairpin", i, j),
            (Gap::None, [(p, q), (j, i)]) =>
                write!(f, "{:<8} ({:>3}, {:>3}), ({:>3}, {:>3})", "Interior", i, j, p, q),
            (Gap::None, [.., (j, i)]) => write!(f, "{:<8} ({:>3}, {:>3}), {}",
                "Multibr.", i, j, fmt_branches(&canonical[..canonical.len() - 1])),
            (Gap::None, []) => unreachable!("A closed loop has a closing pair."),
            (Gap::Open { five, three }, _) if !joint => write!(f, "{:<8} [{:>3}, {:>3}], {}",
                "Exterior", five, three, fmt_branches(&canonical)),
            (Gap::Open { five, three }, _) => write!(f, "{:<8} [{:>3}, {:>3}], {}",
                "Jxterior", five, three, fmt_branches(&canonical)),
            (Gap::Wrap { last, first }, _) => write!(f, "{:<8} [{:>3}, {:>3}], {}",
                "Circular", first, last, fmt_branches(&canonical)),
        }
    }
}

/// The operations of the move sets, for the loops of one strand without
/// pseudoknots, linear or circular: pairs are added and removed within one
/// loop at a time. Joint exterior loops and the loops of a pseudoknot are
/// not covered (their pairs are not in index order).
impl NestedLoop {
    /// The closing pair (i, j) of a closed loop.
    pub fn closing(&self) -> Option<(usize, usize)> {
        let (j, i) = *self.pairs().last().filter(|_| self.gap == Gap::None)?;
        Some((i, j))
    }

    /// All pairs (i, j) with i < j: the closing pair first, then in loop order.
    pub fn index_pairs(&self) -> Vec<(usize, usize)> {
        let (branches, _) = self.branches();
        let mut pairs = Vec::with_capacity(branches.len() + 1);
        pairs.extend(self.closing());
        pairs.extend_from_slice(branches);
        pairs
    }

    /// The first and the last index of the loop: the closing pair, or the
    /// strand ends.
    pub fn span(&self) -> (usize, usize) {
        match self.gap {
            Gap::None => self.closing().unwrap(),
            Gap::Open { five, three } => (five, three),
            Gap::Wrap { last, first } => (first, last),
        }
    }

    /// The stretches of the loop between helix ends, in index order: from
    /// the start of the loop (its closing pair or 5' end) to the first
    /// branch, between branches, and from the last branch to the end.
    pub fn stretches(&self) -> impl Iterator<Item = Stretch> + '_ {
        let (start, end) = self.span();
        let (branches, _) = self.branches();
        // An open loop or a circle starts and ends at unpaired bases.
        let ends = self.gap != Gap::None;
        let n = branches.len();
        let from = std::iter::once(start).chain(branches.iter().map(|&(_, q)| q));
        let to = branches.iter().map(|&(p, _)| p).chain(std::iter::once(end));
        from.zip(to).enumerate().map(move |(k, (a, b))| Stretch {
            ends: (a, b),
            unpaired: if ends && k == 0 { a } else { a + 1 }..if ends && k == n { b + 1 } else { b },
        })
    }

    /// The number of unpaired bases of a hairpin: a closed loop with one
    /// pair, or the loop through the wrap with one pair.
    pub fn hairpin_size(&self) -> Option<usize> {
        match (self.gap, self.pairs()) {
            (Gap::None, [(j, i)]) => Some(j - i - 1),
            (Gap::Wrap { last, first }, [(p, q)]) => Some((last - q) + (p - first)),
            _ => None,
        }
    }

    /// Two stacked pairs: a closed loop, or the loop through the wrap, with
    /// two pairs and no unpaired base. Decided by index, for loops of nested
    /// structures, which pass at most the wrap of their own circle; the
    /// stacks of a pseudoknot are `Pseudoknot::stacks`.
    pub fn is_stack(&self) -> bool {
        match (self.gap, self.pairs()) {
            (Gap::None, [(p, q), (j, i)]) => *p == i + 1 && *j == q + 1,
            (Gap::Wrap { last, first }, [(p, q), (r, s)]) => *p == first && *s == last && *r == q + 1,
            _ => false,
        }
    }

    /// The branches in index order, and the closing pair of a closed loop.
    fn branches(&self) -> (&[Pair], &[Pair]) {
        let pairs = self.pairs();
        match self.gap {
            Gap::None => pairs.split_at(pairs.len() - 1),
            _ => (pairs, &[]),
        }
    }

    /// Add the pair (i, j) of two unpaired bases of the loop: the outer loop
    /// keeps the gap, the inner loop is closed by (i, j). Returns (outer, inner).
    pub fn split_loop(&self, i: usize, j: usize) -> (Self, Self) {
        debug_assert!(i < j, "Split pair (i,j) must satisfy i < j");
        let (branches, closing) = self.branches();
        let lo = branches.partition_point(|&(a, _)| a < i);
        let hi = branches.partition_point(|&(a, _)| a < j);
        debug_assert!(branches[lo..hi].iter().all(|&(_, b)| b < j), "Pair ({i}, {j}) crosses a branch.");
        let outer = Pairs::concat(&[&branches[..lo], &[(i, j)], &branches[hi..], closing]);
        let inner = Pairs::concat(&[&branches[lo..hi], &[(j, i)]]);
        (Self { pairs: outer, gap: self.gap }, Self { pairs: inner, gap: Gap::None })
    }

    /// Remove the closing pair of `inner`, a branch of this loop: its pairs
    /// take the place of that branch.
    pub fn join_loop(&self, inner: &Self) -> Self {
        let (j, i) = inner.closing().map(|(i, j)| (j, i)).expect("The inner loop is closed.");
        let (inside, _) = inner.branches();
        let pairs = self.pairs();
        let k = self.branches().0.partition_point(|&(a, _)| a < i);
        debug_assert_eq!(pairs[k], (i, j), "The inner loop is a branch.");
        let pairs = Pairs::concat(&[&pairs[..k], inside, &pairs[k + 1..]]);
        Self { pairs, gap: self.gap }
    }

    /// The open loop one base longer at its 3' end.
    pub fn extended(&self) -> Self {
        let Gap::Open { five, three } = self.gap else {
            panic!("Only an open loop is extended.");
        };
        Self { pairs: self.pairs.clone(), gap: Gap::Open { five, three: three + 1 } }
    }
}

/// Branches as "[  i,   j], ..." for display.
fn fmt_branches(branches: &[(usize, usize)]) -> String {
    branches.iter()
        .map(|(i, j)| format!("[{:>3}, {:>3}]", i, j))
        .collect::<Vec<_>>()
        .join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NucleotideVec;

    fn lp(pairs: &[(usize, usize)], gap: Gap) -> NestedLoop {
        NestedLoop::new(Pairs::from(pairs.to_vec()), gap)
    }

    #[test]
    fn split_and_join() {
        let open = Gap::Open { five: 0, three: 40 };
        let wrap = Gap::Wrap { last: 40, first: 0 };
        for (l, (i, j), outer, inner) in [
            (lp(&[(20, 2)], Gap::None), (5, 12), lp(&[(5, 12), (20, 2)], Gap::None), lp(&[(12, 5)], Gap::None)),
            (lp(&[(10, 20), (30, 0)], Gap::None), (3, 8),
                lp(&[(3, 8), (10, 20), (30, 0)], Gap::None), lp(&[(8, 3)], Gap::None)),
            (lp(&[(10, 20), (30, 0)], Gap::None), (5, 25),
                lp(&[(5, 25), (30, 0)], Gap::None), lp(&[(10, 20), (25, 5)], Gap::None)),
            (lp(&[(10, 20), (30, 0)], Gap::None), (22, 28),
                lp(&[(10, 20), (22, 28), (30, 0)], Gap::None), lp(&[(28, 22)], Gap::None)),
            (lp(&[(2, 8), (10, 18), (22, 38), (40, 0)], Gap::None), (9, 20),
                lp(&[(2, 8), (9, 20), (22, 38), (40, 0)], Gap::None), lp(&[(10, 18), (20, 9)], Gap::None)),
            (lp(&[(2, 8), (10, 18)], open), (1, 30), lp(&[(1, 30)], open), lp(&[(2, 8), (10, 18), (30, 1)], Gap::None)),
            (lp(&[], open), (5, 12), lp(&[(5, 12)], open), lp(&[(12, 5)], Gap::None)),
            (lp(&[(5, 20)], wrap), (25, 35), lp(&[(5, 20), (25, 35)], wrap), lp(&[(35, 25)], Gap::None)),
            (lp(&[], wrap), (5, 12), lp(&[(5, 12)], wrap), lp(&[(12, 5)], Gap::None)),
        ] {
            assert_eq!(l.split_loop(i, j), (outer.clone(), inner.clone()), "{l:?} + ({i}, {j})");
            assert_eq!(outer.join_loop(&inner), l, "{l:?} - ({i}, {j})");
        }
        // The unpaired bases of a multiloop.
        let multi = lp(&[(2, 8), (9, 20), (22, 38), (40, 0)], Gap::None);
        assert_eq!(multi.index_pairs(), vec![(0, 40), (2, 8), (9, 20), (22, 38)]);
        assert_eq!(multi.stretches().flat_map(|s| s.unpaired).collect::<Vec<_>>(), vec![1, 21, 39]);
    }

    #[test]
    fn move_set_operations() {
        let ends = Gap::Open { five: 0, three: 12 };
        let wrap = Gap::Wrap { last: 12, first: 0 };
        // (loop, closing, index pairs, span, stretches as (ends, unpaired))
        for (l, closing, index_pairs, span, stretches) in [
            (lp(&[(9, 2)], Gap::None), Some((2, 9)), vec![(2, 9)], (2, 9), vec![((2, 9), 3..9)]),
            (lp(&[(4, 7), (9, 2)], Gap::None), Some((2, 9)), vec![(2, 9), (4, 7)], (2, 9),
                vec![((2, 4), 3..4), ((7, 9), 8..9)]),
            (lp(&[(2, 5), (7, 10), (12, 0)], Gap::None), Some((0, 12)), vec![(0, 12), (2, 5), (7, 10)], (0, 12),
                vec![((0, 2), 1..2), ((5, 7), 6..7), ((10, 12), 11..12)]),
            (lp(&[(2, 5), (7, 10)], ends), None, vec![(2, 5), (7, 10)], (0, 12),
                vec![((0, 2), 0..2), ((5, 7), 6..7), ((10, 12), 11..13)]),
            (lp(&[(2, 5), (7, 10)], wrap), None, vec![(2, 5), (7, 10)], (0, 12),
                vec![((0, 2), 0..2), ((5, 7), 6..7), ((10, 12), 11..13)]),
            (lp(&[], ends), None, vec![], (0, 12), vec![((0, 12), 0..13)]),
        ] {
            assert_eq!(l.closing(), closing, "{l:?}");
            assert_eq!(l.index_pairs(), index_pairs, "{l:?}");
            assert_eq!(l.span(), span, "{l:?}");
            let expected: Vec<Stretch> = stretches.into_iter().map(|(ends, unpaired)| Stretch { ends, unpaired }).collect();
            assert_eq!(l.stretches().collect::<Vec<_>>(), expected, "{l:?}");
        }
        // Hairpins, also through the wrap; stacks, also through the wrap.
        assert_eq!(lp(&[(9, 2)], Gap::None).hairpin_size(), Some(6));
        assert_eq!(lp(&[(2, 10)], wrap).hairpin_size(), Some(4));
        assert_eq!(lp(&[(4, 7), (9, 2)], Gap::None).hairpin_size(), None);
        assert!(lp(&[(3, 8), (9, 2)], Gap::None).is_stack());
        assert!(!lp(&[(4, 7), (9, 2)], Gap::None).is_stack());
        assert!(lp(&[(0, 5), (6, 12)], wrap).is_stack());
        assert!(!lp(&[(2, 5), (7, 10)], wrap).is_stack());
        // The exterior loop of a growing strand.
        assert_eq!(lp(&[(2, 5)], ends).extended(), lp(&[(2, 5)], Gap::Open { five: 0, three: 13 }));
    }

    #[test]
    fn segments_of_loops() {
        let seq = NucleotideVec::try_from("GGGAAACCCAGGGAAACCCAAAGGGAAACCCAA").unwrap();
        let ends = Gap::Open { five: 0, three: 32 };
        for (l, ranges) in [
            (lp(&[(6, 2)], Gap::None), vec![2..=6]),
            (lp(&[(2, 6), (7, 1)], Gap::None), vec![1..=2, 6..=7]),
            (lp(&[(2, 8), (10, 18), (22, 28), (30, 0)], Gap::None), vec![0..=2, 8..=10, 18..=22, 28..=30]),
            (lp(&[(2, 8), (10, 18)], ends), vec![0..=2, 8..=10, 18..=32]),
            (lp(&[], ends), vec![0..=32]),
        ] {
            let expected: Vec<&[Base]> = ranges.into_iter().map(|r| &seq[r]).collect();
            let segments: Vec<_> = l.segments(&seq).collect();
            assert_eq!(segments.iter().map(|s| s.as_ref()).collect::<Vec<_>>(), expected);
        }
        // The hairpin through the wrap: 30..=32 and 0..=2.
        let circ = lp(&[(2, 30)], Gap::Wrap { last: 32, first: 0 });
        assert_eq!(circ.segments(&seq).collect::<Vec<_>>(), vec![[&seq[30..=32], &seq[0..=2]].concat()]);
    }
}
