//! Pseudoknots in the ring decomposition. All loops are traced, also those
//! that meet their pairs out of order. The loop graph has a vertex per
//! loop and an edge per pair, between the loops on its two sides.
//! Without pseudoknots it is a forest. A pseudoknot is made of its blocks
//! (biconnected components): the pairs on a cycle, and the loops they join.

use std::fmt;

use ff_structure::MultiPairTable;
use ff_structure::MultiStruct;

use crate::Gap;
use crate::NestedLoop;
use super::nested_loops::Pairs;
use super::rings::Ring;

/// A loop of a pseudoknot, other than an interior loop of a strained stem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PseudoknotLoop {
    pub nested: NestedLoop,
    /// Per segment, true if strained: two visits of the same block enclose
    /// it. The segment before each pair, and for an open loop the segment
    /// after its last pair (the order of `NestedLoop::segments`).
    pub strained: Vec<bool>,
    /// Per pair, true if it ends a strained stem (it lies on a cycle).
    pub stem_ends: Vec<bool>,
}

/// A pseudoknot: the loops of blocks of the loop graph that share
/// loops. Everything else is a regular loop.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Pseudoknot {
    pub loops: Vec<PseudoknotLoop>,
    /// The stacks of the strained stems: two pairs without an unpaired base
    /// between them along the backbone, also through the wraps of circles.
    pub stacks: Vec<NestedLoop>,
    /// The other interior loops of the strained stems: bulges and interior
    /// loops.
    pub interior_loops: Vec<NestedLoop>,
    /// The number of strained stems: chains of pairs on a cycle, joined by
    /// the loops in `stacks` and `interior_loops`.
    pub stems: usize,
    /// The number of knots: blocks, joined where they share strained segments
    /// of a loop (always in a closed loop).
    pub knots: usize,
    /// The number of independent cycles of the loop graph in its blocks:
    /// pairs minus loops plus one, per block.
    pub cycles: usize,
}

impl fmt::Display for Pseudoknot {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:<8} {} loops, {} knots, {} strained stems",
            "Pseudok.", self.loops.len(), self.knots, self.stems)?;
        // With "{:#}", also list its loops.
        if f.alternate() {
            for l in &self.loops {
                write!(f, "\n  {}", l.nested)?;
            }
            for l in self.stacks.iter().chain(&self.interior_loops) {
                write!(f, "\n  {l}")?;
            }
        }
        Ok(())
    }
}

const NONE: usize = usize::MAX;

/// Walk along the backbone from `start`: at a paired base, continue after
/// its partner; at a strand wrap, continue at the first base of the circle.
/// Marks the paired bases it arrives at and the wraps it passes with `id`.
/// The walk ends at a strand end (an open loop) or where it began: at a
/// marked paired base (a closed loop) or at a marked wrap (the loop through
/// the wrap, marked before the walk starts at the first base of the circle).
/// Returns the loop and the number of unpaired bases it passed.
fn walk(mpt: &MultiPairTable, start: usize, id: usize, arrive: &mut [usize]) -> (NestedLoop, usize) {
    let n = mpt.len();
    let mut pairs = Pairs::None;
    let mut unpaired = 0;
    let mut c = start;
    loop {
        // The end of the structure ends its last strand.
        match (c < n).then(|| mpt.get(c)) {
            None | Some(MultiStruct::StrandBreak) =>
                return (NestedLoop::new(pairs, Gap::Open { five: start, three: c - 1 }), unpaired),
            Some(MultiStruct::Unpaired) => {
                unpaired += 1;
                c += 1;
            }
            Some(MultiStruct::StrandWrap) if arrive[c] != NONE =>
                return (NestedLoop::new(pairs, Gap::Wrap { last: c - 1, first: start }), unpaired),
            Some(MultiStruct::StrandWrap) => {
                arrive[c] = id;
                c = first_of_circle(mpt, c);
            }
            Some(MultiStruct::Paired(_)) if arrive[c] != NONE => {
                // Back at the start: the closing pair, the one seen from its
                // largest index, goes last.
                let last = pairs.as_slice().iter().enumerate().max_by_key(|(_, p)| p.0).unwrap().0;
                pairs.rotate_left(last + 1);
                return (NestedLoop::new(pairs, Gap::None), unpaired);
            }
            Some(MultiStruct::Paired(q)) => {
                arrive[c] = id;
                pairs.push((c, q));
                c = q + 1;
            }
        }
    }
}

/// The first base of the circle that ends with the strand wrap `w`: after
/// the strand end before it.
fn first_of_circle(mpt: &MultiPairTable, w: usize) -> usize {
    (0..w).rev()
        .find(|&k| matches!(mpt.get(k), MultiStruct::StrandBreak | MultiStruct::StrandWrap))
        .map_or(0, |e| e + 1)
}

/// The loop graph of a structure: a vertex per loop, an edge per pair.
#[derive(Debug)]
struct LoopGraph {
    /// The loops, each with the unpaired bases along its walk.
    loops: Vec<(NestedLoop, usize)>,
    /// The pairs (i, j), i < j, in the order of i.
    pairs: Vec<(usize, usize)>,
    /// Per pair, the loops on its two sides: the loop that arrives at i and
    /// the loop that arrives at j.
    edges: Vec<(usize, usize)>,
}

/// The loop graph of a structure, with all its loops. The open loops first,
/// each walked from the 5' end of its strand once the strand ends; then the
/// loops through the wraps that no open loop passes; then the closed loops.
fn trace_loops(mpt: &MultiPairTable) -> LoopGraph {
    let n = mpt.len();
    let mut arrive = vec![NONE; n];
    let mut loops = Vec::new();
    let mut circles = Vec::new();
    // The first base of the current strand.
    let mut first = 0;
    for k in 0..n {
        match mpt.get(k) {
            MultiStruct::StrandBreak => loops.push(walk(mpt, first, loops.len(), &mut arrive)),
            MultiStruct::StrandWrap => circles.push((first, k)),
            _ => continue,
        }
        first = k + 1;
    }
    if first < n {
        loops.push(walk(mpt, first, loops.len(), &mut arrive));
    }
    for (first, w) in circles {
        if arrive[w] == NONE {
            arrive[w] = loops.len();
            loops.push(walk(mpt, first, loops.len(), &mut arrive));
        }
    }
    for a in 0..n {
        if matches!(mpt.get(a), MultiStruct::Paired(_)) && arrive[a] == NONE {
            loops.push(walk(mpt, a, loops.len(), &mut arrive));
        }
    }
    let (mut pairs, mut edges) = (Vec::new(), Vec::new());
    for i in 0..n {
        if let MultiStruct::Paired(j) = mpt.get(i) && i < j {
            pairs.push((i, j));
            edges.push((arrive[i], arrive[j]));
        }
    }
    LoopGraph { loops, pairs, edges }
}

/// The biconnected components (Hopcroft & Tarjan) of a multigraph, as a
/// block index for each edge that lies on a cycle, None for bridges. The
/// depth-first search keeps its path on a stack.
fn blocks(n: usize, edges: &[(usize, usize)]) -> Vec<Option<usize>> {
    let mut adjacent = vec![Vec::new(); n];
    for (e, &(u, v)) in edges.iter().enumerate() {
        adjacent[u].push((e, v));
        adjacent[v].push((e, u));
    }
    let mut found = vec![NONE; n];
    let mut low = vec![0; n];
    let mut block = vec![None; edges.len()];
    let mut blocks = 0;
    // The edges of the blocks not yet closed.
    let mut stack = Vec::new();
    // The search path: each vertex, the edge it was reached by, and the
    // index of its next adjacent edge.
    let mut path: Vec<(usize, usize, usize)> = Vec::new();
    let mut time = 0;
    for start in 0..n {
        if found[start] != NONE {
            continue;
        }
        found[start] = time;
        low[start] = time;
        time += 1;
        path.push((start, NONE, 0));
        while let Some(&(u, by, next)) = path.last() {
            let Some(&(e, v)) = adjacent[u].get(next) else {
                // Back from u to its parent w.
                path.pop();
                if let Some(&(w, _, _)) = path.last() {
                    low[w] = low[w].min(low[u]);
                    if low[u] >= found[w] {
                        // w separates the block that contains the edge to u,
                        // unless that edge is a bridge.
                        let at = stack.iter().rposition(|&f| f == by).unwrap();
                        let edges_of_block = stack.split_off(at);
                        if edges_of_block.len() > 1 {
                            edges_of_block.into_iter().for_each(|f| block[f] = Some(blocks));
                            blocks += 1;
                        }
                    }
                }
                continue;
            };
            path.last_mut().unwrap().2 += 1;
            if e == by {
                continue;
            }
            if found[v] == NONE {
                stack.push(e);
                found[v] = time;
                low[v] = time;
                time += 1;
                path.push((v, e, 0));
            } else if v == u {
                // A self-loop is a block of its own (each one is seen twice).
                if block[e].is_none() {
                    block[e] = Some(blocks);
                    blocks += 1;
                }
            } else if found[v] < found[u] {
                stack.push(e);
                low[u] = low[u].min(found[v]);
            }
        }
    }
    block
}

fn root(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn join(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (root(parent, a), root(parent, b));
    parent[ra] = rb;
}

/// The first and the last visit of a block to a loop.
struct Span { first: usize, last: usize, block: usize }

/// The span of each block of a loop, in the order of first visits.
/// `blocks` holds the block of each pair, in loop order.
fn spans(blocks: &[Option<usize>]) -> Vec<Span> {
    let mut spans: Vec<Span> = Vec::new();
    for (t, b) in blocks.iter().enumerate() {
        let Some(b) = *b else { continue };
        match spans.iter_mut().find(|s| s.block == b) {
            Some(s) => s.last = t,
            None => spans.push(Span { first: t, last: t, block: b }),
        }
    }
    spans
}

/// Per segment of a loop (see `PseudoknotLoop::strained`), true if two
/// visits of the same block enclose it. In a closed loop (or the loop
/// through a wrap), every segment is enclosed; in an open loop, those from
/// the first visit of a block to its last. `blocks` holds the block of each
/// pair, in loop order.
fn strained_segments(gap: Gap, blocks: &[Option<usize>]) -> Vec<bool> {
    let k = blocks.len();
    if !matches!(gap, Gap::Open { .. }) {
        return vec![true; k];
    }
    let mut strained = vec![false; k + 1];
    for s in spans(blocks) {
        (s.first + 1..=s.last).for_each(|t| strained[t] = true);
    }
    strained
}

/// The blocks of a loop that share a strained segment, as pairs to join
/// into one knot. In a closed loop, all of its blocks; in an open loop, the
/// blocks of overlapping spans.
fn shared_blocks(gap: Gap, blocks: &[Option<usize>]) -> Vec<(usize, usize)> {
    if !matches!(gap, Gap::Open { .. }) {
        let mut present = blocks.iter().flatten();
        let Some(&a) = present.next() else { return Vec::new() };
        return present.filter(|&&b| b != a).map(|&b| (a, b)).collect();
    }
    // In order of first visits; `reach` is the end of the overlapping spans
    // so far, and one of their blocks.
    let mut shared = Vec::new();
    let mut reach: Option<(usize, usize)> = None;
    for s in spans(blocks) {
        reach = match reach {
            Some((end, b)) if s.first < end => {
                shared.push((b, s.block));
                Some((end.max(s.last), s.block))
            }
            _ => Some((s.last, s.block)),
        };
    }
    shared
}

/// Emit the rings of a structure with crossing pairs: the loops that lie on
/// no cycle of the loop graph as `Ring::Loop`, then the pseudoknots.
pub(super) fn trace_pseudoknots<F: FnMut(&Ring)>(mpt: &MultiPairTable, f: &mut F) {
    let LoopGraph { loops, pairs, edges } = trace_loops(mpt);
    let block = blocks(loops.len(), &edges);
    // The block of each pair, at both of its bases.
    let mut block_at = vec![None; mpt.len()];
    for (&(i, j), &b) in pairs.iter().zip(&block) {
        (block_at[i], block_at[j]) = (b, b);
    }
    // Loops that share a block belong to one pseudoknot.
    let mut group: Vec<usize> = (0..loops.len()).collect();
    for (&(u, v), b) in edges.iter().zip(&block) {
        if b.is_some() {
            join(&mut group, u, v);
        }
    }
    let mut knots: Vec<usize> = (0..block.iter().flatten().max().map_or(0, |b| b + 1)).collect();

    /// A pseudoknot while its loops are assigned.
    #[derive(Default)]
    struct Part {
        pk: Pseudoknot,
        /// The blocks of its loops, once per loop.
        blocks: Vec<usize>,
        /// The visits of its loops to pairs on a cycle: twice per pair.
        visits: usize,
    }
    let mut parts: Vec<Part> = Vec::new();
    let mut part_of = vec![NONE; loops.len()];
    for (v, (l, unpaired)) in loops.into_iter().enumerate() {
        let visits: Vec<Option<usize>> = l.pairs().iter().map(|&(a, _)| block_at[a]).collect();
        if visits.iter().all(Option::is_none) {
            f(&Ring::Loop(l));
            continue;
        }
        let g = root(&mut group, v);
        if part_of[g] == NONE {
            part_of[g] = parts.len();
            parts.push(Part::default());
        }
        let part = &mut parts[part_of[g]];
        let mut blocks: Vec<usize> = visits.iter().flatten().copied().collect();
        part.visits += blocks.len();
        blocks.sort_unstable();
        blocks.dedup();
        part.blocks.extend(blocks);
        // Two pairs on a cycle, in a closed loop or the loop through the wrap:
        // a stack without unpaired bases, else an interior loop.
        let closed = matches!(l.gap(), Gap::None | Gap::Wrap { .. });
        if closed && visits.len() == 2 && visits.iter().all(Option::is_some) {
            match unpaired {
                0 => part.pk.stacks.push(l),
                _ => part.pk.interior_loops.push(l),
            }
        } else {
            shared_blocks(l.gap(), &visits).into_iter().for_each(|(a, b)| join(&mut knots, a, b));
            let strained = strained_segments(l.gap(), &visits);
            let stem_ends = visits.iter().map(Option::is_some).collect();
            part.pk.loops.push(PseudoknotLoop { nested: l, strained, stem_ends });
        }
    }
    for Part { mut pk, mut blocks, visits } in parts {
        // Per block: pairs minus loops plus one; the loops of a block are
        // the entries in `blocks` before the dedup.
        let pairs = visits / 2;
        let loops_of_blocks = blocks.len();
        blocks.sort_unstable();
        blocks.dedup();
        pk.cycles = pairs + blocks.len() - loops_of_blocks;
        // Pairs on a cycle, minus the links between pairs of one stem.
        pk.stems = pairs - pk.stacks.len() - pk.interior_loops.len();
        let mut roots: Vec<usize> = blocks.into_iter().map(|b| root(&mut knots, b)).collect();
        roots.sort_unstable();
        roots.dedup();
        pk.knots = roots.len();
        f(&Ring::Pseudoknot(pk));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RingDecomposition;

    #[test]
    fn loop_graphs() {
        // An H-type: the open loop and the stacks of the two stems, one
        // cycle per stem through the open loop.
        let g = trace_loops(&MultiPairTable::try_from("((([[[)))...]]]").unwrap());
        assert_eq!(g.loops.len(), 5);
        assert_eq!(g.loops[0].0.gap(), Gap::Open { five: 0, three: 14 });
        assert_eq!(g.loops.iter().map(|l| l.1).collect::<Vec<_>>(), [3, 0, 0, 0, 0]);
        assert_eq!(g.pairs, [(0, 8), (1, 7), (2, 6), (3, 14), (4, 13), (5, 12)]);
        assert_eq!(g.edges, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]);
        // The same, circular: the open loop becomes the loop through the wrap.
        let c = trace_loops(&MultiPairTable::try_from("((([[[)))...]]]^").unwrap());
        assert_eq!(c.loops[0].0.gap(), Gap::Wrap { last: 14, first: 0 });
        assert_eq!(c.edges, g.edges);
        // A circle joined to a linear strand: two edges between the same two
        // loops, a cycle.
        let g = trace_loops(&MultiPairTable::try_from("((..+..))^").unwrap());
        assert_eq!(g.loops.len(), 2);
        assert_eq!(g.loops[0].0.pairs(), [(0, 8), (7, 1)]);
        assert_eq!(g.edges, [(0, 1), (1, 0)]);
    }

    #[test]
    fn strained_and_shared() {
        let (a, b, c, d) = (Some(0), Some(1), Some(2), Some(3));
        let open = Gap::Open { five: 0, three: 99 };
        // Crossing spans (an H-type): strained between them, one knot.
        assert_eq!(strained_segments(open, &[a, b, a, b]), [false, true, true, true, false]);
        assert_eq!(shared_blocks(open, &[a, b, a, b]), [(0, 1)]);
        // Two H-types side by side: the linker between them is free.
        let two = [a, b, a, b, c, d, c, d];
        assert_eq!(strained_segments(open, &two), [false, true, true, true, false, true, true, true, false]);
        assert_eq!(shared_blocks(open, &two), [(0, 1), (2, 3)]);
        // A span inside another: one knot. Spans that only touch: two.
        assert_eq!(shared_blocks(open, &[a, b, b, a]), [(0, 1)]);
        assert_eq!(strained_segments(open, &[a, a, b, b]), [false, true, false, true, false]);
        assert!(shared_blocks(open, &[a, a, b, b]).is_empty());
        // A branch between two visits of a block is enclosed.
        assert_eq!(strained_segments(open, &[a, None, a]), [false, true, true, false]);
        // A closed loop: every segment, and all of its blocks share them.
        assert_eq!(strained_segments(Gap::None, &[a, None, b]), [true, true, true]);
        assert_eq!(shared_blocks(Gap::None, &[a, None, b, a]), [(0, 1)]);
        assert!(shared_blocks(Gap::None, &[a, a]).is_empty());
    }

    #[test]
    fn blocks_of_multigraphs() {
        // A triangle with a pendant edge, and a pair of parallel edges.
        let edges = [(0, 1), (1, 2), (2, 0), (2, 3), (4, 5), (5, 4)];
        let block = blocks(6, &edges);
        assert_eq!(block[0], block[1]);
        assert_eq!(block[1], block[2]);
        assert!(block[0].is_some());
        assert_eq!(block[3], None);
        assert_eq!(block[4], block[5]);
        assert!(block[4].is_some() && block[4] != block[0]);
        // A self-loop is on a cycle, a single edge is not.
        assert!(blocks(2, &[(0, 0), (0, 1)])[0].is_some());
        assert_eq!(blocks(2, &[(0, 0), (0, 1)])[1], None);
        // A long cycle, deeper than a recursive search could go.
        let n = 200_000;
        let ring: Vec<(usize, usize)> = (0..n).map(|u| (u, (u + 1) % n)).collect();
        assert!(blocks(n, &ring).iter().all(|&b| b == Some(0)));
    }

    fn pseudoknots(db: &str) -> Vec<Pseudoknot> {
        MultiPairTable::try_from(db).unwrap().rings().into_iter()
            .filter_map(|r| match r { Ring::Pseudoknot(pk) => Some(pk), _ => None })
            .collect()
    }

    /// Per loop of a pseudoknot: its gap and which segments are strained.
    fn strained(pk: &Pseudoknot) -> Vec<(Gap, Vec<bool>)> {
        pk.loops.iter().map(|l| (l.nested.gap(), l.strained.clone())).collect()
    }

    #[test]
    fn h_type() {
        // One open loop: the tails lead to the strand ends.
        let [pk] = &pseudoknots("((([[[)))...]]]")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems, pk.stacks.len(), pk.interior_loops.len(), pk.cycles), (1, 2, 4, 0, 2));
        assert_eq!(strained(pk), vec![
            (Gap::Open { five: 0, three: 14 }, vec![false, true, true, true, false])]);
        assert!(pk.loops[0].stem_ends.iter().all(|&e| e));
    }

    #[test]
    fn kissing_hairpins() {
        // One strand: the linker between the stems is strained, the tails
        // are not; the pseudoloop is closed.
        let [pk] = &pseudoknots("((((..[[[[..))))((((..]]]]..))))")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems, pk.cycles), (1, 3, 2));
        assert_eq!(strained(pk), vec![
            (Gap::Open { five: 0, three: 31 }, vec![false, true, false]),
            (Gap::None, vec![true; 4])]);

        // Enclosed by a pair: the enclosing loop is strained, the closing
        // pair is a branch of it.
        let [pk] = &pseudoknots("(.((((..[[[[..))))((((..]]]]..)))).)")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems), (1, 3));
        let enclosing = pk.loops.iter().find(|l| l.nested.closing() == Some((0, 35))).unwrap();
        assert_eq!(enclosing.strained, vec![true; 3]);
        assert_eq!(enclosing.stem_ends, vec![true, true, false]);

        // A nick in the linker: the stems are bridges, only the kiss and the
        // pseudoloop remain.
        let [pk] = &pseudoknots("((((..[[[[..))))+((((..]]]]..))))")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems, pk.cycles), (1, 1, 1));
        assert_eq!(strained(pk), vec![(Gap::None, vec![true; 4])]);
        assert_eq!(pk.loops[0].stem_ends, vec![true, false, true, false]);
    }

    #[test]
    fn stems_through_the_wrap() {
        // A circular H-type, cut inside a stem: the stack through the wrap
        // continues the stem, and the loop through the wrap is strained.
        let [pk] = &pseudoknots("((([[[)))...]]]^")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems, pk.stacks.len(), pk.interior_loops.len()), (1, 2, 4, 0));
        let [pk] = &pseudoknots("(([[[))(...]]])^")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems, pk.stacks.len(), pk.interior_loops.len()), (1, 2, 4, 0));
        assert!(pk.stacks.iter().any(|l| matches!(l.gap(), Gap::Wrap { .. })));
        assert!(pk.loops.iter().all(|l| l.strained.iter().all(|&s| s)));
        // A bulge in a strained stem is an interior loop.
        let [pk] = &pseudoknots("((.(([[[))))...]]]")[..] else { panic!() };
        assert_eq!((pk.stems, pk.stacks.len(), pk.interior_loops.len()), (2, 4, 1));
    }

    #[test]
    fn stacks_through_the_wraps_of_two_circles() {
        // A stem of three pairs between two circles, drawn with consecutive
        // indices, and cut inside the stem in both circles: the backbone
        // passes both wraps between two of its pairs. Both are two stacks.
        for db in ["..........(........+.........[[[).....^...]]]............^",
                   "([..............((^))...............)^..........]........"] {
            let [pk] = &pseudoknots(db)[..] else { panic!("{db}") };
            assert_eq!((pk.stems, pk.stacks.len(), pk.interior_loops.len()), (2, 2, 0), "{db}");
        }
    }

    #[test]
    fn knots_of_one_loop() {
        // Two H-types next to each other share the open loop, not its linker.
        let [pk] = &pseudoknots("((([[[)))...]]]..((([[[)))...]]]")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems), (2, 4));
        let [(_, segments)] = &strained(pk)[..] else { panic!() };
        assert_eq!(segments, &[false, true, true, true, false, true, true, true, false]);
        // An H-type in the pseudoloop of another one: one knot.
        let [pk] = &pseudoknots("((([[[)))..<<<{{{>>>}}}..]]]")[..] else { panic!() };
        assert_eq!((pk.knots, pk.stems), (1, 4));
    }
}
