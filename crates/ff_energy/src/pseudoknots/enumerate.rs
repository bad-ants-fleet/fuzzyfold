//! Band finding algorithm by Condon et al.
//! main enumeration file

use ff_structure::{ExtendedDotBracketVec, PairTable, StructureError};
use std::collections::HashSet;

use crate::pseudoknots::{
    ClosedRegion, ClosingDescriptor, LocationStatus, Loop, LoopType, PseudoloopContext, RegionTree,
    build_closed_regions_tree, closing_descriptor, closing_pairs, collect_bands,
    extended_dot_bracket_to_pair_table, interior_loop_type, is_pseudo, location_status,
    nested_pairs,
};

/// Enumerates the SPAN_BAND loops for `region`: the interior loops (or
/// multiloops, if other children sit in the gaps) between consecutive
/// rungs of each band found by `collect_bands`.
pub fn enumerate_band_spanning_loops(
    tree: &RegionTree,
    region: &ClosedRegion,
    pt: &PairTable,
    loops: &mut Vec<Loop>,
) {
    let bands = collect_bands(tree, region, pt);

    for chain in &bands {
        for idx in 0..chain.len().saturating_sub(1) {
            let outer_left = chain[idx];
            let inner_left = chain[idx + 1];
            let outer_right = pt[outer_left].expect("band rung must be paired") as usize;
            let inner_right = pt[inner_left].expect("band rung must be paired") as usize;

            let closing_bp = (outer_left, outer_right);
            let inner_bp = (inner_left, inner_right);

            let mut nest = nested_pairs(tree, &region.children, pt, outer_left, inner_left);
            nest.extend(nested_pairs(
                tree,
                &region.children,
                pt,
                inner_right,
                outer_right,
            ));

            if nest.is_empty() {
                // No nested closing pairs in the gap: Stack / Bulge / Interior.
                let (ltype, n5, n3) = interior_loop_type(closing_bp, inner_bp);
                loops.push(
                    Loop::new(ltype, LocationStatus::SpanBand)
                        .with_closing(ClosingDescriptor::Single(closing_bp))
                        .with_inner(inner_bp)
                        .with_unpaired(n5, n3),
                );
            } else {
                // A nested closing pair sits between the two consecutive band rungs →
                // the segment is a SpanBand Multiloop. The next band rung (inner_bp) is
                // placed first in children, followed by any incidental nested pairs.
                // From the energy model's perspective all children are equivalent branches;
                // bp (DP feature 8) is charged once per branch regardless of which branch
                // is the band continuation.
                //
                // The unpaired counts feed the cp term (DP feature 9), which charges only
                // bases lying directly in this multiloop. Bases inside a branch's own loop
                // are excluded — that branch's loop energy already accounts for them.
                let n_unpaired_5p = count_unpaired_outside_children(
                    tree, &region.children, pt, outer_left + 1, inner_left,
                );
                let n_unpaired_3p = count_unpaired_outside_children(
                    tree, &region.children, pt, inner_right + 1, outer_right,
                );
                let mut children = vec![ClosingDescriptor::Single(inner_bp)];
                children.extend(nest.into_iter().map(ClosingDescriptor::Single));

                loops.push(
                    Loop::new(LoopType::Multiloop, LocationStatus::SpanBand)
                        .with_closing(ClosingDescriptor::Single(closing_bp))
                        .with_children(children)
                        .with_unpaired(n_unpaired_5p, n_unpaired_3p),
                );
            }
        }
    }
}

// Main loop enumeration function

/// Builds the list of `Loop`s for a structure, given its closed-regions tree.
/// Mirrors the Python `enumerate_loops(root, pt)`, with the `region.is_root`
/// branch replaced by a final pass over `tree.root_children`.
pub fn enumerate_loops(tree: &RegionTree, pt: &PairTable) -> Vec<Loop> {
    let mut loops = Vec::new();

    for &idx in &tree.root_children {
        visit(tree, idx, pt, &mut loops, PseudoloopContext::External);
    }

    // External loop: in Python this was the `region.is_root` branch of
    // `visit`, using `root.children`. Here `root_children` plays that role.
    let child_bps: Vec<ClosingDescriptor> = tree
        .root_children
        .iter()
        .map(|&idx| closing_descriptor(&tree.nodes[idx], pt))
        .collect();

    loops.push(Loop::new(LoopType::External, LocationStatus::Standard).with_children(child_bps));

    loops
}

/// Post-order visit of `idx` and its descendants, appending one `Loop` per
/// region (plus any SPAN_BAND loops from `enumerate_band_spanning_loops`).
///
/// `pk_context` is the context this node sees — i.e. what type of loop
/// encloses this node. Passed top-down so `Pseudoloop` loops can record the
/// correct initiation-penalty context.
fn visit(
    tree: &RegionTree,
    idx: usize,
    pt: &PairTable,
    loops: &mut Vec<Loop>,
    pk_context: PseudoloopContext,
) {
    let region = &tree.nodes[idx];
    let n_children = region.children.len();
    let pseudo = is_pseudo(region, pt);

    // Determine the context that children of THIS node will inherit.
    let is_multiloop = !pseudo
        && (n_children > 1 || (n_children == 1 && is_pseudo(&tree.nodes[region.children[0]], pt)));
    let child_context = if pseudo {
        PseudoloopContext::Pseudoloop
    } else if is_multiloop {
        PseudoloopContext::Multiloop
    } else {
        pk_context
    };

    for &child_idx in &region.children {
        visit(tree, child_idx, pt, loops, child_context);
    }

    let closing = closing_pairs(region, pt);
    let loc = location_status(region, region.parent.map(|p| &tree.nodes[p]), pt);

    // Multiloops arise in three places in this function:
    //   1. Here, pseudo branch: a pseudoloop's children list may contain multiple
    //      band-tip pairs and nested regions — but the enclosing loop is Pseudoloop,
    //      not Multiloop; the multiloop accounting is handled via n_pup/n_bands/n_nested.
    //   2. Here, else branch (line ~237): a non-pseudo region with ≥2 children, or
    //      1 child that is itself a pseudoloop → LoopType::Multiloop.
    //   3. In enumerate_band_spanning_loops: a gap between two consecutive band rungs
    //      that contains at least one nested closing pair → LoopType::Multiloop,
    //      LocationStatus::SpanBand. Both the next band rung and any incidental nested
    //      pairs become children (the next rung is just another branch from this loop's
    //      perspective; bp is charged once per branch regardless).
    if pseudo {
        let closing_set: HashSet<(usize, usize)> = closing.iter().copied().collect();
        let mut child_pairs: Vec<ClosingDescriptor> = Vec::new();

        let bands: Vec<Vec<usize>> = collect_bands(tree, region, pt);

        for chain in &bands {
            let outer_left_arm = chain[0];
            let inner_left_arm = *chain.last().unwrap();
            let outer_bp = (
                outer_left_arm,
                pt[outer_left_arm].expect("band rung must be paired") as usize,
            );
            let inner_bp = (
                inner_left_arm,
                pt[inner_left_arm].expect("band rung must be paired") as usize,
            );

            if !closing_set.contains(&outer_bp) {
                child_pairs.push(ClosingDescriptor::Single(outer_bp));
            }
            child_pairs.push(ClosingDescriptor::Single(inner_bp));
        }

        for &child_idx in &region.children {
            child_pairs.push(closing_descriptor(&tree.nodes[child_idx], pt));
        }

        child_pairs.sort_by_key(|cd| match cd {
            ClosingDescriptor::Single((i, _)) => *i,
            ClosingDescriptor::Double((i, _), _) => *i,
        });

        // Compute n_pup: free linker nucleotides across all 2n-1 inter-arm gaps.
        //
        // Each band contributes two arm segments (5' arm and 3' arm). Sorting all
        // 2n arm segments by start position and summing the 2n-1 gaps between
        // consecutive segments gives the total span of single-stranded linker
        // nucleotides that belong to the pseudoloop itself. Bases inside nested
        // child regions are subtracted by count_pup_linkers, so each unpaired base
        // is counted in exactly one loop (the Rastegari partition for the pup term).
        //
        // No special case for n=2 — the formula is identical for any band count.
        // Band count is not density (density.rs): a density-2 structure may have
        // any number of bands. What it does bound is the Dirks-Pierce *algorithm*
        // ("precisely two spanning regions", DP 2003); the energy model is not so
        // limited (Andronescu 2010), and this crate only evaluates.
        if bands.len() > 2 {
            // Only shown with RUST_LOG=warn and a logger installed.
            log::warn!(
                "pseudoloop with {} bands. The Dirks-Pierce algorithm emits only \
                 2-band pseudoloops, and pb was fit on a predominantly 2-band \
                 corpus, so this value is weakly constrained. (Band count is not \
                 density: a density-2 structure may have any number of bands.)",
                bands.len()
            );
        }
        let mut arm_segs: Vec<(usize, usize)> = bands.iter().flat_map(|chain| {
            let outer    = chain[0];
            let inner    = *chain.last().unwrap();
            let inner_3p = pt[inner].expect("band rung must be paired") as usize;
            let outer_3p = pt[outer].expect("band rung must be paired") as usize;
            [(outer, inner), (inner_3p, outer_3p)]
        }).collect();
        arm_segs.sort_by_key(|&(s, _)| s);
        let n_pup: usize = arm_segs.windows(2)
            .map(|w| count_pup_linkers(tree, &region.children, pt, w[0].1 + 1, w[1].0))
            .sum();

        loops.push(
            Loop::new(LoopType::Pseudoloop, loc)
                .with_closing(ClosingDescriptor::Double(closing[0], closing[1]))
                .with_children(child_pairs)
                .with_pup(n_pup)
                .with_bands(bands.len())
                .with_nested(n_children)
                .with_pk_context(pk_context),
        );
    } else if n_children == 0 {
        let (ci, cj) = closing[0];
        loops.push(
            Loop::new(LoopType::Hairpin, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_unpaired(cj - ci - 1, 0),
        );
    } else if n_children == 1 && !is_pseudo(&tree.nodes[region.children[0]], pt) {
        let inner_pair = closing_pairs(&tree.nodes[region.children[0]], pt)[0];
        let (ltype, n5, n3) = interior_loop_type(closing[0], inner_pair);
        loops.push(
            Loop::new(ltype, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_inner(inner_pair)
                .with_unpaired(n5, n3),
        );
    } else {
        // Non-pseudo region with ≥2 children, or 1 child that is itself a pseudoloop.
        // Classification matches NearestNeighborLoop::classify: branches.len() >= 2
        // → Multibranch (here LoopType::Multiloop). A single pseudoloop child also
        // lands here because the pseudoloop counts as a branch of the enclosing loop.
        let child_bps: Vec<ClosingDescriptor> = region
            .children
            .iter()
            .map(|&c| closing_descriptor(&tree.nodes[c], pt))
            .collect();
        loops.push(
            Loop::new(LoopType::Multiloop, loc)
                .with_closing(ClosingDescriptor::Single(closing[0]))
                .with_children(child_bps),
        );
    }

    enumerate_band_spanning_loops(tree, region, pt, loops);
}

/// Count unpaired positions in `[start, end)` that lie directly in the enclosing
/// loop — i.e. excluding any position inside a nested child region.
///
/// This is the `cp` term of a SpanBand multiloop (DP feature 9). A base inside a
/// branch's hairpin or interior loop belongs to that branch, not to the multiloop,
/// and is already charged by the branch's own loop energy; counting it here too
/// would charge it twice. Mirrors the child-span exclusion that
/// [`count_pup_linkers`] performs for the pseudoloop `pup` term.
fn count_unpaired_outside_children(
    tree: &RegionTree,
    children: &[usize],
    pt: &PairTable,
    start: usize,
    end: usize,
) -> usize {
    if start >= end {
        return 0;
    }
    let child_spans: Vec<(usize, usize)> = children
        .iter()
        .filter_map(|&ci| {
            let pairs = closing_pairs(&tree.nodes[ci], pt);
            let cb = pairs.iter().map(|&(i, _)| i).min()?;
            let ce = pairs.iter().map(|&(_, j)| j).max()?;
            Some((cb, ce))
        })
        .collect();
    (start..end)
        .filter(|&p| {
            pt[p].is_none() && !child_spans.iter().any(|&(cb, ce)| cb <= p && p <= ce)
        })
        .count()
}

/// Count free linker positions in the pseudoloop gap `[start, end)`.
///
/// Mirrors HotKnots' `pseudoEnergyDP`: takes the raw positional span of the
/// gap and subtracts the full span of every nested child region that overlaps
/// with it.  This matches the physical intent of the DP09 `pup` parameter —
/// only genuinely single-stranded linker nucleotides between structured
/// elements are counted, not bases that belong to a child hairpin or stem.
fn count_pup_linkers(
    tree: &RegionTree,
    children: &[usize],
    pt: &PairTable,
    start: usize,
    end: usize,
) -> usize {
    if start >= end {
        return 0;
    }
    let raw = end - start;
    let child_overlap: usize = children
        .iter()
        .filter_map(|&ci| {
            let pairs = closing_pairs(&tree.nodes[ci], pt);
            if pairs.is_empty() {
                return None;
            }
            let cb = pairs.iter().map(|&(i, _)| i).min().unwrap();
            let ce = pairs.iter().map(|&(_, j)| j).max().unwrap();
            // Intersection of child span [cb, ce] with gap [start, end-1].
            let lo = cb.max(start);
            let hi = ce.min(end - 1);
            if lo <= hi { Some(hi - lo + 1) } else { None }
        })
        .sum();
    raw - child_overlap
}

/// Builds the list of `Loop`s for a structure directly from its `PairTable`,
/// skipping the string round-trip.
///
/// Prefer this over [`parse_structure`] whenever the pair table is already
/// available (e.g. during beam search): it avoids the O(P²)
/// `pair_table_to_dot_bracket` step and the redundant re-parse back to a
/// `PairTable` that `parse_structure` performs internally.
pub fn parse_loops_from_pt(pt: &PairTable) -> Vec<Loop> {
    let tree = build_closed_regions_tree(pt);
    enumerate_loops(&tree, pt)
}

/// Main entry point for the whole enumeration process.
/// Takes a dot-bracket string, builds the pair table and closed-regions tree,
/// and returns the list of `Loop`s.
pub fn parse_structure(s: &str) -> Result<Vec<Loop>, StructureError> {
    let edbv = ExtendedDotBracketVec::try_from(s)?;
    let pt = extended_dot_bracket_to_pair_table(&edbv)?;
    Ok(parse_loops_from_pt(&pt))
}

// TESTS

#[cfg(test)]
mod band_spanning_loop_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_h_type_pseudoknot_band_spanning_loops() {
        // ((([[[)))...]]] -> bands [[0,1,2],[3,4,5]], no children
        // -> 4 stacked-pair (Stack) loops, all SPAN_BAND
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let tree = build_closed_regions_tree(&pt);
        let region = &tree.nodes[tree.root_children[0]];

        let mut loops = Vec::new();
        enumerate_band_spanning_loops(&tree, region, &pt, &mut loops);

        assert_eq!(loops.len(), 4);
        for lp in &loops {
            assert_eq!(lp.loop_type, LoopType::Stack);
            assert_eq!(lp.location, LocationStatus::SpanBand);
            assert_eq!((lp.unpaired_5p, lp.unpaired_3p), (0, 0));
        }

        assert_eq!(loops[0].closing, Some(ClosingDescriptor::Single((0, 8))));
        assert_eq!(loops[0].inner, Some((1, 7)));

        assert_eq!(loops[1].closing, Some(ClosingDescriptor::Single((1, 7))));
        assert_eq!(loops[1].inner, Some((2, 6)));

        assert_eq!(loops[2].closing, Some(ClosingDescriptor::Single((3, 14))));
        assert_eq!(loops[2].inner, Some((4, 13)));

        assert_eq!(loops[3].closing, Some(ClosingDescriptor::Single((4, 13))));
        assert_eq!(loops[3].inner, Some((5, 12)));
    }

    #[test]
    fn test_loop_display() {
        let lp = Loop::new(LoopType::Stack, LocationStatus::SpanBand)
            .with_closing(ClosingDescriptor::Single((0, 8)))
            .with_inner((1, 7));

        assert_eq!(
            lp.to_string(),
            "Loop(Stack, SpanBand, closing=(0, 8), inner=(1, 7))"
        );
    }
}

#[cfg(test)]
mod enumerate_loops_tests {
    use super::*;

    #[test]
    fn test_simple_hairpin() {
        // (((...)))
        let loops = parse_structure("(((...)))").unwrap();
        assert_eq!(loops.len(), 4);

        assert_eq!(
            loops[0],
            Loop::new(LoopType::Hairpin, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((2, 6)))
                .with_unpaired(3, 0)
        );

        assert_eq!(
            loops[1],
            Loop::new(LoopType::Stack, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((1, 7)))
                .with_inner((2, 6))
        );

        assert_eq!(
            loops[2],
            Loop::new(LoopType::Stack, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Single((0, 8)))
                .with_inner((1, 7))
        );

        assert_eq!(
            loops[3],
            Loop::new(LoopType::External, LocationStatus::Standard)
                .with_children(vec![ClosingDescriptor::Single((0, 8))])
        );
    }

    #[test]
    fn test_h_type_pseudoknot() {
        // ((([[[)))...]]]
        let loops = parse_structure("((([[[)))...]]]").unwrap();
        assert_eq!(loops.len(), 6);

        assert_eq!(
            loops[0],
            Loop::new(LoopType::Pseudoloop, LocationStatus::Standard)
                .with_closing(ClosingDescriptor::Double((0, 8), (3, 14)))
                .with_children(vec![
                    ClosingDescriptor::Single((2, 6)),
                    ClosingDescriptor::Single((5, 12)),
                ])
                .with_pup(3)
                .with_bands(2)
                .with_nested(0)
                .with_pk_context(PseudoloopContext::External)
        );

        for lp in &loops[1..5] {
            assert_eq!(lp.loop_type, LoopType::Stack);
            assert_eq!(lp.location, LocationStatus::SpanBand);
        }

        assert_eq!(
            loops[5],
            Loop::new(LoopType::External, LocationStatus::Standard)
                .with_children(vec![ClosingDescriptor::Double((0, 8), (3, 14))])
        );
    }
}

/// Structural invariants of the loop enumeration.
///
/// These pin properties that follow from Rastegari & Condon's decomposition
/// rather than from any one structure: that the unpaired bases partition
/// across loops, that every `LocationStatus` and `PseudoloopContext` variant
/// is reachable and under what conditions, and how band segments are
/// classified. Several were written after the corresponding behaviour turned
/// out to be un-exercised — `OutBand` in particular cannot arise unless the
/// gap between nested border pairs is wide enough to hold a closed region.
#[cfg(test)]
mod enumeration_invariants {
    use super::*;

    fn find(loops: &[Loop], lt: LoopType, closing: (usize, usize)) -> &Loop {
        let cd = ClosingDescriptor::Single(closing);
        loops
            .iter()
            .find(|l| l.loop_type == lt && l.closing == Some(cd))
            .unwrap_or_else(|| panic!("no {lt:?} loop with closing {closing:?}"))
    }

    fn pseudoloop(loops: &[Loop]) -> &Loop {
        loops
            .iter()
            .find(|l| l.loop_type == LoopType::Pseudoloop)
            .expect("structure has no pseudoloop")
    }

    // ── Rastegari & Condon partition invariant ───────────────────────────

    /// Unpaired positions in `[lo, hi]` that lie outside every child span.
    fn free_in(pt: &PairTable, children: &[ClosingDescriptor], lo: usize, hi: usize) -> usize {
        let spans: Vec<(usize, usize)> = children
            .iter()
            .map(|cd| match cd {
                ClosingDescriptor::Single((a, b)) => (*a, *b),
                ClosingDescriptor::Double((a1, _), (_, b2)) => (*a1, *b2),
            })
            .collect();
        (lo..=hi)
            .filter(|&q| pt[q].is_none() && !spans.iter().any(|&(a, b)| a <= q && q <= b))
            .count()
    }

    /// Rastegari & Condon: *every unpaired base is in exactly one loop*.
    ///
    /// Each loop type is asked how many unpaired bases it owns, using the same
    /// accounting the energy model bills for. Summing over all loops must
    /// reproduce the structure's unpaired count exactly: a double-count makes
    /// the sum too high, an unowned base makes it too low.
    fn assert_unpaired_partition(dbn: &str) {
        let edbv = ExtendedDotBracketVec::try_from(dbn).unwrap();
        let pt = extended_dot_bracket_to_pair_table(&edbv).unwrap();
        let loops = parse_structure(dbn).unwrap();

        let total = (0..pt.len()).filter(|&q| pt[q].is_none()).count();

        let owned: usize = loops
            .iter()
            .map(|lp| match lp.loop_type {
                LoopType::Hairpin
                | LoopType::Stack
                | LoopType::Bulge
                | LoopType::Interior => lp.unpaired_5p + lp.unpaired_3p,
                LoopType::Pseudoloop => lp.n_pup,
                LoopType::Multiloop => {
                    if lp.location == LocationStatus::SpanBand {
                        lp.unpaired_5p + lp.unpaired_3p
                    } else {
                        let (i, j) = match lp.closing.unwrap() {
                            ClosingDescriptor::Single(p) => p,
                            ClosingDescriptor::Double(p1, _) => p1,
                        };
                        free_in(&pt, &lp.children, i + 1, j - 1)
                    }
                }
                LoopType::External => free_in(&pt, &lp.children, 0, pt.len() - 1),
            })
            .sum();

        assert_eq!(
            owned, total,
            "unpaired partition violated for {dbn}: loops own {owned}, structure has {total}"
        );
    }

    #[test]
    fn unpaired_partition_holds() {
        for dbn in [
            // nested
            "(((...)))",
            "(.(...)..(...)....)",
            "((..((..))..))",
            "(((...(((...)))...(((...)))...)))",
            // H-type and friends
            "((([[[)))...]]]",
            "(((...[[[...)))...]]]",
            "((([[[))).(...).]]]",
            "((([[[.))).(...).(((.]]])))",
            "(.(..).([[[))...]]]",
            // higher density
            "((((...[[[[...))))....((((...]]]]...))))",
            "(((...[[[...{{{...)))...]]]...}}}",
            "(((...[[[...{{{...)))...]]]...<<<...}}}...>>>",
            // PK in multiloop / PK in PK
            "(((.((((....)))).((((...[[[[...))))...]]]].)))",
            "(((...[[[...)))...{{{...<<<...}}}...>>>...]]]",
            // a 5-pair band chain crossing a single-pair band, branch in one segment
            "(((...(.(..)..(...............[...................)........))))............]",
            // multiloop spanning a band with branches in BOTH gaps
            "(.(...)....(((...[[[..))).(...)..]]]...)",
        ] {
            assert_unpaired_partition(dbn);
        }
    }

    // ── LocationStatus: when each variant arises ─────────────────────────

    /// A parent that is not pseudoknotted (pi(0) = 8) makes every loop in the
    /// structure `Standard`.
    #[test]
    fn standard_when_parent_not_pseudoknotted() {
        let loops = parse_structure("(((...)))").unwrap();
        assert!(
            loops.iter().all(|l| l.location == LocationStatus::Standard),
            "non-pseudoknotted structure must be entirely Standard"
        );
    }

    /// The H-type pseudoknot `((([[[)))...]]]` yields one SpanBand loop per
    /// consecutive rung pair per band — 4 in total, all `Stack` with
    /// n5 = n3 = 0.
    #[test]
    fn spanband_loop_per_consecutive_rung_pair() {
        let loops = parse_structure("((([[[)))...]]]").unwrap();

        let span: Vec<_> = loops
            .iter()
            .filter(|l| l.location == LocationStatus::SpanBand)
            .collect();
        assert_eq!(span.len(), 4, "4 SpanBand loops: 2 per band");
        for lp in &span {
            assert_eq!(lp.loop_type, LoopType::Stack);
            assert_eq!((lp.unpaired_5p, lp.unpaired_3p), (0, 0));
        }

        let pk = pseudoloop(&loops);
        assert_eq!(pk.n_bands, 2);
        assert_eq!(pk.n_pup, 3, "the 3 unpaired linker bases at 9..11");
    }

    /// When the border pairs CROSS, I1 union I2 covers the whole region and every
    /// child is `InBand`.
    ///
    /// The child must be a real closed region: in `((([[[)))...]]]` positions
    /// 9..11 are unpaired, so nothing there can be classified at all. The gap
    /// is widened here to hold an actual hairpin (10,14).
    #[test]
    fn inband_when_border_pairs_cross() {
        let loops = parse_structure("((([[[))).(...).]]]").unwrap();

        let hp = find(&loops, LoopType::Hairpin, (10, 14));
        assert_eq!(
            hp.location,
            LocationStatus::InBand,
            "child of a pseudoknotted parent with crossing borders is InBand"
        );

        // Borders cross: pi(0) = 8 >= 3 = pi(18).
        let pk = pseudoloop(&loops);
        assert_eq!(pk.closing, Some(ClosingDescriptor::Double((0, 8), (3, 18))));
        assert_eq!(pk.n_bands, 2);
    }

    /// When the border pairs are NESTED rather than crossing, the gap between
    /// them is `OutBand`.
    ///
    /// This variant is easy to believe is covered when it is not. In
    /// `((([[[.)))..(((.]]])))` the gap is only positions 10..11 — too small to
    /// hold a closed region — so that structure produces no OutBand loop at
    /// all. The gap here is widened to fit a hairpin (11,15), which is the
    /// smallest structure that actually reaches the OutBand branch.
    #[test]
    fn outband_needs_a_gap_wide_enough_for_a_child() {
        let loops = parse_structure("((([[[.))).(...).(((.]]])))").unwrap();

        let hp = find(&loops, LoopType::Hairpin, (11, 15));
        assert_eq!(
            hp.location,
            LocationStatus::OutBand,
            "child in the gap between nested border pairs is OutBand"
        );

        // Borders are nested: pi(0) = 9 < 17 = pi(26), leaving the gap (9,17).
        let pk = pseudoloop(&loops);
        assert_eq!(
            pk.closing,
            Some(ClosingDescriptor::Double((0, 9), (17, 26)))
        );
    }

    /// A hairpin closing strictly inside the 5' gap of a band segment makes that
    /// segment a `Multiloop` while its `LocationStatus` stays `SpanBand` —
    /// LoopType and LocationStatus are independent axes.
    ///
    /// The 5' gap spans positions 1..6, but only 1 and 6 belong to this
    /// multiloop: 3 and 4 lie inside the hairpin (2,5) and are charged by that
    /// hairpin's own loop energy instead. DP feature 9 therefore counts 2.
    #[test]
    fn branch_in_gap_promotes_spanband_to_multiloop() {
        let loops = parse_structure("(.(..).([[[))...]]]").unwrap();

        let ml = find(&loops, LoopType::Multiloop, (0, 12));
        assert_eq!(ml.location, LocationStatus::SpanBand);
        assert_eq!(
            ml.children,
            vec![
                ClosingDescriptor::Single((7, 11)), // band continuation
                ClosingDescriptor::Single((2, 5)),  // incidental branch
            ],
            "next band rung first, then nested branches"
        );
        assert_eq!((ml.unpaired_5p, ml.unpaired_3p), (2, 0));

        let hp = find(&loops, LoopType::Hairpin, (2, 5));
        assert_eq!(hp.location, LocationStatus::InBand);
    }

    // ── PseudoloopContext: all three are reachable ───────────────────────

    /// The pseudoloop initiation penalty is context dependent; all three contexts
    /// are reachable. The strings below are minimal valid realisations of each.
    #[test]
    fn every_pseudoloop_context_is_reachable() {
        // Badge 1: a pseudoknot at top level -> init_external.
        let ext = parse_structure("(((...[[[...)))...]]]").unwrap();
        assert_eq!(
            pseudoloop(&ext).pk_context,
            Some(PseudoloopContext::External)
        );

        // Badge 2: a pseudoknot as one branch of a multiloop -> init_multiloop.
        let ml = parse_structure("(((.((((....)))).((((...[[[[...))))...]]]].)))").unwrap();
        assert_eq!(
            pseudoloop(&ml).pk_context,
            Some(PseudoloopContext::Multiloop)
        );

        // Badge 3: a pseudoknot nested inside another pseudoknot's linker
        // -> init_pseudoloop. The inner PK uses {} and <> so its bands stay
        // distinct from the outer PK's () and [].
        let pk = parse_structure("(((...[[[...)))...{{{...<<<...}}}...>>>...]]]").unwrap();
        let inner = pk
            .iter()
            .find(|l| {
                l.loop_type == LoopType::Pseudoloop
                    && l.closing == Some(ClosingDescriptor::Double((18, 32), (24, 38)))
            })
            .expect("inner pseudoloop");
        assert_eq!(inner.pk_context, Some(PseudoloopContext::Pseudoloop));

        // The outer pseudoloop counts the inner one as a nested region.
        let outer = pk
            .iter()
            .find(|l| {
                l.loop_type == LoopType::Pseudoloop
                    && l.closing == Some(ClosingDescriptor::Double((0, 14), (6, 44)))
            })
            .expect("outer pseudoloop");
        assert_eq!(outer.pk_context, Some(PseudoloopContext::External));
        assert_eq!(outer.n_nested, 1);
    }

    // ── Band counts across pseudoknot topologies ─────────────────────────

    /// Density (the most pairs crossing any one point) is NOT the same as
    /// `n_bands` (the number of crossing helices). H-type is the only topology
    /// where the two coincide, which makes the distinction easy to miss.
    #[test]
    fn band_counts_across_pseudoknot_topologies() {
        // H-type — density 2, 2 bands.
        let h = parse_structure("(((...[[[...)))...]]]").unwrap();
        assert_eq!(pseudoloop(&h).n_bands, 2);

        // K-type kissing hairpin — density 2, but 3 bands: the
        // bridging stem crosses both hairpin stems, which do not cross
        // each other.
        let k = parse_structure("((((...[[[[...))))....((((...]]]]...))))").unwrap();
        assert_eq!(pseudoloop(&k).n_bands, 3);

        // L-type 3-knot — density 3, 3 mutually crossing helices.
        let l = parse_structure("(((...[[[...{{{...)))...]]]...}}}").unwrap();
        assert_eq!(pseudoloop(&l).n_bands, 3);

        // M-type 4-knot — density 3 (NOT 4), 4 bands: A,B,C mutually
        // cross and D crosses only C.
        let m = parse_structure("(((...[[[...{{{...)))...]]]...<<<...}}}...>>>").unwrap();
        assert_eq!(pseudoloop(&m).n_bands, 4);
    }

    // ── A plain multiloop, for contrast with the SpanBand case ───────────

    /// A plain multiloop has the same anatomy as a SpanBand one — one closing
    /// pair, two branches, unpaired bases directly in the loop — and differs
    /// only in which parameter set the energy model applies. Bases inside the
    /// branches' own hairpin loops are not counted here; they belong to those
    /// hairpins.
    #[test]
    fn plain_multiloop_has_no_spanband_loops() {
        let loops = parse_structure("(.(...)..(...)....)").unwrap();

        let ml = find(&loops, LoopType::Multiloop, (0, 18));
        assert_eq!(
            ml.location,
            LocationStatus::Standard,
            "no band involved, so this is a plain multiloop"
        );
        assert_eq!(
            ml.children,
            vec![
                ClosingDescriptor::Single((2, 6)),
                ClosingDescriptor::Single((9, 13)),
            ]
        );

        // Unpaired bases belonging directly to the loop: 1, 7, 8 and 14..17.
        for closing in [(2, 6), (9, 13)] {
            let hp = find(&loops, LoopType::Hairpin, closing);
            assert_eq!(hp.location, LocationStatus::Standard);
            assert_eq!((hp.unpaired_5p, hp.unpaired_3p), (3, 0));
        }
        assert!(
            loops.iter().all(|l| l.location != LocationStatus::SpanBand),
            "nothing here spans a band"
        );
    }

    // ── Band segment classification ──────────────────────────────────────

    /// One structure exercising every band-segment classification at once.
    /// Band A is a 5-pair chain (0,62) -> (1,61) -> (2,60) -> (6,59) -> (14,50)
    /// crossing the single-pair Band B (30,75), with a hairpin (8,11) nested
    /// in the last band segment. The four segment gaps demonstrate every
    /// SpanBand classification:
    ///
    ///   (0,62)->(1,61)   empty gap                  -> Stack     (feature 10)
    ///   (1,61)->(2,60)   empty gap                  -> Stack     (feature 10)
    ///   (2,60)->(6,59)   unpaired only, one-sided   -> Bulge     (feature 11)
    ///   (6,59)->(14,50)  contains a branch          -> Multiloop (features 7,8,9)
    ///
    /// The nested hairpin (8,11) is InBand (feature 6).
    ///
    /// The third segment is a `Bulge` rather than `Interior`: it is the
    /// one-sided special case, the gap being empty on the 3' side. Both take
    /// the same `e_intp` scaling in the energy model, so the DP feature is the
    /// same either way.
    #[test]
    fn band_segments_classify_by_gap_contents() {
        let s = "(((...(.(..)..(...............[...................)........))))............]";
        let loops = parse_structure(s).unwrap();

        let pk = pseudoloop(&loops);
        assert_eq!(pk.n_bands, 2, "Band A (5-pair chain) and Band B (1 pair)");

        // Feature 10: two empty segments -> stacked pairs spanning a band.
        for closing in [(0, 62), (1, 61)] {
            let lp = find(&loops, LoopType::Stack, closing);
            assert_eq!(lp.location, LocationStatus::SpanBand);
        }

        // Feature 11: unpaired bases only on the 5' side -> Bulge.
        let bulge = find(&loops, LoopType::Bulge, (2, 60));
        assert_eq!(bulge.location, LocationStatus::SpanBand);
        assert_eq!((bulge.unpaired_5p, bulge.unpaired_3p), (3, 0));

        // Features 7/8/9: the segment holding a branch -> SpanBand Multiloop.
        let ml = find(&loops, LoopType::Multiloop, (6, 59));
        assert_eq!(ml.location, LocationStatus::SpanBand);
        assert_eq!(
            ml.children,
            vec![
                ClosingDescriptor::Single((14, 50)), // band continuation
                ClosingDescriptor::Single((8, 11)),  // the nested hairpin
            ]
        );
        // The cp term (DP feature 9) counts 7, 12 and 13 on the 5' side.
        // Bases 9 and 10 lie inside the hairpin (8,11) and belong to its loop.
        assert_eq!((ml.unpaired_5p, ml.unpaired_3p), (3, 8));

        // Feature 6: the nested substructure sits inside a band.
        let hp = find(&loops, LoopType::Hairpin, (8, 11));
        assert_eq!(hp.location, LocationStatus::InBand);
    }
}
