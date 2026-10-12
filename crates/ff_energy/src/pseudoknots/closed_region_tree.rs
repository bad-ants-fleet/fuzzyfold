//! Build the tree of closed regions from a pair table (Algorithm 1,
//! Rastegari & Condon), adapted to 0-based indexing and an arena
//! representation with no explicit root node.

use ff_structure::PairTable;
use std::collections::{HashMap, HashSet};
use std::fmt;

/// A node in the closed-regions tree. `i` and `j` are 0-based endpoints.
///
/// `parent == None` means this region is top-level — i.e. its Python
/// counterpart's `parent` was the sentinel root `(-1, n)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClosedRegion {
    pub i: usize,
    pub j: usize,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
}

impl ClosedRegion {
    fn new(i: usize, j: usize) -> Self {
        ClosedRegion {
            i,
            j,
            parent: None,
            children: Vec::new(),
        }
    }
}

/// Arena of `ClosedRegion`s. There is no node for the root; `n` and
/// `root_children` replace the Python root's `j` and `children`.
#[derive(Debug, Clone)]
pub struct RegionTree {
    pub nodes: Vec<ClosedRegion>,
    /// Length of the structure (the implicit root's right endpoint).
    pub n: usize,
    /// Direct children of the implicit root (arena indices, sorted by `i`).
    pub root_children: Vec<usize>,
}

/// Re-parents `region_idx` as a root child, claiming any currently-pooled
/// regions nested inside it (those with `i > region.i`).
///
/// Mirrors `add_to_tree` from the Python reference implementation, where
/// `root.children` acted as the pool; here `root_children` is that pool.
fn add_to_tree(nodes: &mut [ClosedRegion], root_children: &mut Vec<usize>, region_idx: usize) {
    let region_i = nodes[region_idx].i;

    let pooled = std::mem::take(root_children);
    let mut new_children = Vec::new();
    let mut remaining = Vec::new();

    for child_idx in pooled {
        if nodes[child_idx].i > region_i {
            new_children.push(child_idx);
        } else {
            remaining.push(child_idx);
        }
    }

    new_children.sort_by_key(|&c| nodes[c].i);
    for &child_idx in &new_children {
        nodes[child_idx].parent = Some(region_idx);
    }

    nodes[region_idx].children = new_children;
    // nodes[region_idx].parent stays None — it's top-level.

    remaining.push(region_idx);
    *root_children = remaining;
}

/// Builds the tree of closed regions from a `PairTable` (Algorithm 1,
/// Rastegari & Condon, adapted to 0-indexing, no sentinel root node).
pub fn build_closed_regions_tree(pt: &PairTable) -> RegionTree {
    let n = pt.len();

    let mut nodes: Vec<ClosedRegion> = Vec::new();
    let mut root_children: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();

    for lam in 0..n {
        if let Some(b) = pt[lam] {
            let b = b as usize;

            if lam < b {
                // Case 1: opening bracket
                let idx = nodes.len();
                nodes.push(ClosedRegion::new(lam, b));
                stack.push(idx);
            } else if b < lam {
                // Case 2: closing bracket — merge crossing (pseudoknotted) regions
                let mut e = lam;
                while let Some(&top) = stack.last() {
                    if nodes[top].i > b {
                        e = e.max(nodes[top].j);
                        stack.pop();
                    } else {
                        break;
                    }
                }
                if let Some(&top) = stack.last() {
                    nodes[top].j = nodes[top].j.max(e);
                }
            }
        }

        // Case 3: does lam close the region on top of the stack?
        if let Some(&top) = stack.last()
            && lam == nodes[top].j
        {
            let region_idx = stack.pop().unwrap();
            add_to_tree(&mut nodes, &mut root_children, region_idx);
        }
    }

    RegionTree {
        nodes,
        n,
        root_children,
    }
}

// Display implementations for debugging and visualization
impl fmt::Display for ClosedRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ClosedRegion[{},{}]", self.i, self.j)
    }
}

impl fmt::Display for RegionTree {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "root (n={})", self.n)?;

        let last = self.root_children.len().saturating_sub(1);
        for (idx, &child) in self.root_children.iter().enumerate() {
            self.fmt_node(f, child, "", idx == last)?;
        }
        Ok(())
    }
}

impl RegionTree {
    fn fmt_node(
        &self,
        f: &mut fmt::Formatter<'_>,
        idx: usize,
        prefix: &str,
        is_last: bool,
    ) -> fmt::Result {
        let region = &self.nodes[idx];
        let connector = if is_last { "└── " } else { "├── " };
        writeln!(f, "{prefix}{connector}{region}")?;

        let new_prefix = format!("{prefix}{}", if is_last { "    " } else { "│   " });
        let last = region.children.len().saturating_sub(1);
        for (cidx, &child) in region.children.iter().enumerate() {
            self.fmt_node(f, child, &new_prefix, cidx == last)?;
        }
        Ok(())
    }
}

// for python export and visualizations
// ── Step-by-step instrumented build ──────────────────────────────────────────

#[derive(Debug, Clone)]
pub enum BuildEvent {
    Unpaired,
    Opening,
    Closing { merged: Vec<(usize, usize)> },
}

#[derive(Debug, Clone)]
pub struct BuildStep {
    pub lam: usize,
    pub event: BuildEvent,
    /// Arena indices of open regions (not yet completed).
    pub stack: Vec<usize>,
    /// Arena indices of completed root children (same pool as `RegionTree::root_children`).
    pub root_children: Vec<usize>,
    /// Snapshot of the full arena after this step.
    pub nodes: Vec<ClosedRegion>,
    /// Arena index of the region completed this step (Case 3), if any.
    pub completed: Option<usize>,
}

/// Same algorithm as `build_closed_regions_tree`, but emits one `BuildStep`
/// per `lam` iteration so callers can observe the tree growing incrementally.
pub fn build_closed_regions_tree_steps(pt: &PairTable) -> Vec<BuildStep> {
    let n = pt.len();
    let mut nodes: Vec<ClosedRegion> = Vec::new();
    let mut root_children: Vec<usize> = Vec::new();
    let mut stack: Vec<usize> = Vec::new();
    let mut steps: Vec<BuildStep> = Vec::new();

    for lam in 0..n {
        let mut event = BuildEvent::Unpaired;
        let mut completed = None;

        if let Some(b) = pt[lam] {
            let b = b as usize;

            if lam < b {
                // Case 1: opening bracket
                let idx = nodes.len();
                nodes.push(ClosedRegion::new(lam, b));
                stack.push(idx);
                event = BuildEvent::Opening;
            } else if b < lam {
                // Case 2: closing bracket — merge crossing regions
                let mut merged = Vec::new();
                let mut e = lam;
                while let Some(&top) = stack.last() {
                    if nodes[top].i > b {
                        e = e.max(nodes[top].j);
                        merged.push((nodes[top].i, nodes[top].j));
                        stack.pop();
                    } else {
                        break;
                    }
                }
                if let Some(&top) = stack.last() {
                    nodes[top].j = nodes[top].j.max(e);
                }
                event = BuildEvent::Closing { merged };
            }
        }

        // Case 3: does lam close the region on top of the stack?
        if let Some(&top) = stack.last()
            && lam == nodes[top].j
        {
            let region_idx = stack.pop().unwrap();
            add_to_tree(&mut nodes, &mut root_children, region_idx);
            completed = Some(region_idx);
        }

        steps.push(BuildStep {
            lam,
            event,
            stack: stack.clone(),
            root_children: root_children.clone(),
            nodes: nodes.clone(),
            completed,
        });
    }

    steps
}

// Helper fucntions

/// A region is pseudoknotted if its border positions `i` and `j` don't pair
/// with each other (i.e. `pt[i] != j`).
pub fn is_pseudo(region: &ClosedRegion, pt: &PairTable) -> bool {
    pt[region.i].map(|b| b as usize) != Some(region.j)
}

/// The closing pair(s) of a region: 1 pair if non-pseudoknotted, 2 pairs
/// (the crossing borders) if pseudoknotted.
pub fn closing_pairs(region: &ClosedRegion, pt: &PairTable) -> Vec<(usize, usize)> {
    let (i, j) = (region.i, region.j);

    if !is_pseudo(region, pt) {
        vec![(i, j)]
    } else {
        // Both i and j are guaranteed paired by the tree-building invariant
        // (they were pushed/extended only via Some(...) pair-table entries).
        let pi = pt[i].expect("region.i must be paired") as usize;
        let pj = pt[j].expect("region.j must be paired") as usize;
        vec![(i, pi), (pj, j)]
    }
}

/// The closing pair(s) of a region, shaped for use as `Loop` children:
/// either a single pair, or a pair of pairs for a pseudoknotted region's
/// two crossing borders.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClosingDescriptor {
    Single((usize, usize)),
    Double((usize, usize), (usize, usize)),
}

impl fmt::Display for ClosingDescriptor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClosingDescriptor::Single((i, j)) => write!(f, "({i}, {j})"),
            ClosingDescriptor::Double((i, j), (k, l)) => write!(f, "(({i}, {j}), ({k}, {l}))"),
        }
    }
}

pub fn closing_descriptor(region: &ClosedRegion, pt: &PairTable) -> ClosingDescriptor {
    let pairs = closing_pairs(region, pt);
    match pairs.len() {
        1 => ClosingDescriptor::Single(pairs[0]),
        2 => ClosingDescriptor::Double(pairs[0], pairs[1]),
        _ => unreachable!("closing_pairs always returns 1 or 2 pairs"),
    }
}

/// Where a region sits relative to its parent's pairing structure.
///
/// `SpanBand` is not produced by `location_status` — per the original
/// implementation, it's assigned later in the main enumeration function,
/// since it depends on the band-finding algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocationStatus {
    Standard,
    SpanBand,
    InBand,
    OutBand,
}

/// Determines whether `region` is in-band, out-band, or standard, based on
/// its position relative to `parent` and `parent`'s pairing status.
///
/// `parent == None` corresponds to the Python `parent is None or
/// parent.is_root` check — in our arena, a top-level region's `parent`
/// field is `None` (its Python counterpart's parent was the sentinel root).
pub fn location_status(
    region: &ClosedRegion,
    parent: Option<&ClosedRegion>,
    pt: &PairTable,
) -> LocationStatus {
    let parent = match parent {
        None => return LocationStatus::Standard,
        Some(p) => p,
    };

    if !is_pseudo(parent, pt) {
        return LocationStatus::Standard;
    }

    // Parent is pseudoknotted: check whether region.i falls in one of the
    // two bands defined by the parent's (non-crossing-with-each-other)
    // closing-pair endpoints.
    let (pi, pj) = (parent.i, parent.j);
    let p_bi = pt[pi].expect("parent.i must be paired") as usize;
    let p_bj = pt[pj].expect("parent.j must be paired") as usize;
    let ri = region.i;

    // [pi, p_bi] and [p_bj, pj] cover [pi, pj] completely only when the
    // border pairs cross (p_bi >= p_bj) — the standard 2-band H-type case,
    // where every child is IN_BAND. When p_bi < p_bj (e.g. kissing loops,
    // where the borders are nested rather than crossing), the gap
    // (p_bi, p_bj) is OUT_BAND.
    if (pi..=p_bi).contains(&ri) || (p_bj..=pj).contains(&ri) {
        LocationStatus::InBand
    } else {
        LocationStatus::OutBand
    }
}

/// Builds `BL` for `region`, then walks it left-to-right using the §3.2
/// stacking relation. Returns a list of chains; each chain is a band's
/// left-arm positions ordered outer → inner.
///
/// Returns an empty list if `region` is not pseudoknotted.
pub fn collect_bands(tree: &RegionTree, region: &ClosedRegion, pt: &PairTable) -> Vec<Vec<usize>> {
    if !is_pseudo(region, pt) {
        return Vec::new();
    }

    // ----- Step 1: construct BL --------------------------------------------
    // BL = paired positions in [region.i, region.j] that are neither inside
    // a nested closed region nor a closing pair of one.
    let child_ranges: Vec<(usize, usize)> = region
        .children
        .iter()
        .map(|&idx| (tree.nodes[idx].i, tree.nodes[idx].j))
        .collect();

    let mut child_closing: HashSet<(usize, usize)> = HashSet::new();
    for &idx in &region.children {
        child_closing.extend(closing_pairs(&tree.nodes[idx], pt));
    }

    let in_child = |k: usize| child_ranges.iter().any(|&(ci, cj)| ci <= k && k <= cj);

    // Unlike Python's `BL: list[int]` (annotated but not enforced to exclude
    // None), this is `Vec<usize>` — the type itself rules out "None" entries.
    let mut bl: Vec<usize> = Vec::new();
    for k in region.i..=region.j {
        let bk = match pt[k] {
            Some(b) => b as usize,
            None => continue,
        };
        if in_child(k) {
            continue;
        }
        let pair = if k < bk { (k, bk) } else { (bk, k) };
        if child_closing.contains(&pair) {
            continue;
        }
        bl.push(k);
    }

    if bl.is_empty() {
        return Vec::new();
    }

    // O(1) Prev/Next via index.
    let pos: HashMap<usize, usize> = bl.iter().enumerate().map(|(idx, &k)| (k, idx)).collect();

    // ----- Step 2: walk BL, partitioning into bands ------------------------
    let mut bands: Vec<Vec<usize>> = Vec::new();
    let mut i_idx = 0;
    while i_idx < bl.len() {
        let bi = bl[i_idx];
        let pt_bi = pt[bi].expect("BL entries are paired") as usize;

        if bi > pt_bi {
            // right-arm endpoint — skip
            i_idx += 1;
            continue;
        }

        // (bi, pt[bi]) is the outer closing pair of this band.
        let mut chain = vec![bi];
        let mut bi_p = bi;
        let mut bj_p = pt_bi;

        // While Next(bi', BL) == bp(Prev(bj', BL)), step inward.
        loop {
            let i_p = pos[&bi_p];
            let j_p = pos[&bj_p];
            if i_p + 1 >= bl.len() || j_p == 0 {
                break;
            }
            let next_bi = bl[i_p + 1];
            let prev_bj = bl[j_p - 1];
            let pt_prev_bj = pt[prev_bj].expect("BL entries are paired") as usize;

            if next_bi == pt_prev_bj {
                bi_p = next_bi;
                bj_p = prev_bj;
                chain.push(bi_p);
            } else {
                break;
            }
        }

        bands.push(chain);
        i_idx = pos[&bi_p] + 1;
    }

    bands
}

/// Closing pairs of children whose interval lies strictly inside `(left, right)`.
pub fn nested_pairs(
    tree: &RegionTree,
    children: &[usize],
    pt: &PairTable,
    left: usize,
    right: usize,
) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    for &idx in children {
        let c = &tree.nodes[idx];
        if left < c.i && c.j < right {
            result.extend(closing_pairs(c, pt));
        }
    }
    result
}

#[cfg(test)]
mod build_closed_region_tree_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_simple_hairpin() {
        // (((...)))  ->  root_children=[(0,8)] -> (1,7) -> (2,6)
        let pt = PairTable::try_from("(((...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        assert_eq!(tree.n, 9);
        assert_eq!(tree.root_children.len(), 1);

        let r0 = &tree.nodes[tree.root_children[0]];
        assert_eq!((r0.i, r0.j), (0, 8));
        assert_eq!(r0.parent, None);
        assert_eq!(r0.children.len(), 1);

        let r1 = &tree.nodes[r0.children[0]];
        assert_eq!((r1.i, r1.j), (1, 7));
        assert_eq!(r1.children.len(), 1);

        let r2 = &tree.nodes[r1.children[0]];
        assert_eq!((r2.i, r2.j), (2, 6));
        assert!(r2.children.is_empty());
    }

    #[test]
    fn test_h_type_pseudoknot_collapses_to_single_region() {
        // ((([[[)))...]]]  ->  root_children=[(0,14)], pseudoknotted, no children
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let tree = build_closed_regions_tree(&pt);

        assert_eq!(tree.root_children.len(), 1);
        let r0 = &tree.nodes[tree.root_children[0]];
        assert_eq!((r0.i, r0.j), (0, 14));
        assert!(r0.children.is_empty());
        assert_ne!(pt[0 as usize].map(|v| v as usize), Some(r0.j)); // pseudoknotted
    }

    #[test]
    fn test_two_sibling_hairpins() {
        let pt = PairTable::try_from("((..))((..))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        assert_eq!(tree.root_children.len(), 2);
        let r0 = &tree.nodes[tree.root_children[0]];
        let r1 = &tree.nodes[tree.root_children[1]];
        assert_eq!((r0.i, r0.j), (0, 5));
        assert_eq!((r1.i, r1.j), (6, 11));
    }

    // ── add_to_tree unit tests ────────────────────────────────────────────────

    // Test 1: empty pool — the new region becomes the only root child with no
    // children of its own.
    #[test]
    fn test_add_to_tree_empty_pool() {
        let mut nodes = vec![ClosedRegion::new(0, 8)]; // arena idx 0
        let mut root_children: Vec<usize> = Vec::new();

        add_to_tree(&mut nodes, &mut root_children, 0);

        assert_eq!(root_children, vec![0]);
        assert!(nodes[0].children.is_empty());
        assert_eq!(nodes[0].parent, None);
    }

    // Test 2: pool contains a region that is nested inside the new one
    // (its `i` > new region's `i`) — it must be claimed as a child.
    #[test]
    fn test_add_to_tree_claims_nested_region() {
        // arena: idx 0 = [2,6] (inner, already finished), idx 1 = [0,8] (outer)
        let mut nodes = vec![
            ClosedRegion::new(2, 6), // idx 0 — inner, currently in pool
            ClosedRegion::new(0, 8), // idx 1 — outer, being added now
        ];
        let mut root_children = vec![0usize]; // inner is parked as root child

        add_to_tree(&mut nodes, &mut root_children, 1); // add outer

        // outer becomes the only root child
        assert_eq!(root_children, vec![1]);
        // outer claims inner as its child
        assert_eq!(nodes[1].children, vec![0]);
        assert_eq!(nodes[1].parent, None);
        // inner's parent is updated to outer
        assert_eq!(nodes[0].parent, Some(1));
    }

    // Test 3: pool contains one sibling (i < new region's i) and one nested
    // region (i > new region's i). Only the nested one is claimed; the sibling
    // stays in root_children alongside the new region.
    #[test]
    fn test_add_to_tree_keeps_sibling_unclaimed() {
        // arena: idx 0 = [0,5] (sibling), idx 1 = [8,12] (nested), idx 2 = [6,14] (new)
        let mut nodes = vec![
            ClosedRegion::new(0, 5),  // idx 0 — sibling, i=0 < 6 → stays
            ClosedRegion::new(8, 12), // idx 1 — nested,  i=8 > 6 → claimed
            ClosedRegion::new(6, 14), // idx 2 — being added now
        ];
        let mut root_children = vec![0usize, 1usize]; // both parked as root children

        add_to_tree(&mut nodes, &mut root_children, 2);

        // sibling (0) stays; new region (2) is appended after it
        assert_eq!(root_children, vec![0, 2]);
        // new region claims only the nested one
        assert_eq!(nodes[2].children, vec![1]);
        assert_eq!(nodes[2].parent, None);
        assert_eq!(nodes[1].parent, Some(2));
        // sibling is untouched
        assert_eq!(nodes[0].parent, None);
        assert!(nodes[0].children.is_empty());
    }
}

#[cfg(test)]
mod display_closed_region_tree_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_display_simple_hairpin() {
        let pt = PairTable::try_from("(((...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let expected = "\
root (n=9)
└── ClosedRegion[0,8]
    └── ClosedRegion[1,7]
        └── ClosedRegion[2,6]
";
        assert_eq!(tree.to_string(), expected);
    }

    #[test]
    fn test_display_two_siblings() {
        let pt = PairTable::try_from("((..))((..))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let expected = "\
root (n=12)
├── ClosedRegion[0,5]
│   └── ClosedRegion[1,4]
└── ClosedRegion[6,11]
    └── ClosedRegion[7,10]
";
        assert_eq!(tree.to_string(), expected);
    }
}
#[test]
fn test_display_turnip_yellow_rna() {
    let pt = PairTable::try_from("(((((.......)))))(((....[[[[[)))...]]]]]....").unwrap();
    let tree = build_closed_regions_tree(&pt);

    let expected = "\
root (n=44)
├── ClosedRegion[0,16]
│   └── ClosedRegion[1,15]
│       └── ClosedRegion[2,14]
│           └── ClosedRegion[3,13]
│               └── ClosedRegion[4,12]
└── ClosedRegion[17,39]
";
    assert_eq!(tree.to_string(), expected);
}

#[cfg(test)]
mod closing_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_non_pseudo_region() {
        // (((...)))  -> innermost region (2,6), pt[2] == 6
        let pt = PairTable::try_from("(((...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        // walk down to the innermost region
        let r0 = tree.root_children[0];
        let r1 = tree.nodes[r0].children[0];
        let r2 = tree.nodes[r1].children[0];
        let region = &tree.nodes[r2];

        assert_eq!((region.i, region.j), (2, 6));
        assert!(!is_pseudo(region, &pt));
        assert_eq!(closing_pairs(region, &pt), vec![(2, 6)]);
        assert_eq!(
            closing_descriptor(region, &pt),
            ClosingDescriptor::Single((2, 6))
        );
    }

    #[test]
    fn test_pseudo_region() {
        // ((([[[)))...]]]  -> single top-level region (0,14), pt[0] == 8 != 14
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let r0 = tree.root_children[0];
        let region = &tree.nodes[r0];

        assert_eq!((region.i, region.j), (0, 14));
        assert!(is_pseudo(region, &pt));
        assert_eq!(closing_pairs(region, &pt), vec![(0, 8), (3, 14)]);
        assert_eq!(
            closing_descriptor(region, &pt),
            ClosingDescriptor::Double((0, 8), (3, 14))
        );
    }
}

#[cfg(test)]
mod location_status_tests {
    use super::*;

    fn region(i: usize, j: usize, parent: Option<usize>) -> ClosedRegion {
        ClosedRegion {
            i,
            j,
            parent,
            children: Vec::new(),
        }
    }

    #[test]
    fn test_root_children_is_standard() {
        let region = region(0, 8, None);
        let pt = PairTable::new(9);
        assert_eq!(
            location_status(&region, None, &pt),
            LocationStatus::Standard
        );
    }

    #[test]
    fn test_non_pseudo_parent_is_standard() {
        // parent (1,12) with pt[1] == 12 -> not pseudoknotted, regardless of child
        let mut pt = PairTable::new(14);
        pt[1 as usize] = Some(12);
        pt[12 as usize] = Some(1);

        let parent = region(1, 12, None);
        let child = region(4, 9, Some(0));

        assert_eq!(
            location_status(&child, Some(&parent), &pt),
            LocationStatus::Standard
        );
    }

    #[test]
    fn test_in_band_crossing_borders() {
        // parent (0,10), pt[0]=6, pt[10]=3 -> crossing borders (p_bi >= p_bj),
        // so every child is in-band
        let mut pt = PairTable::new(11);
        pt[0 as usize] = Some(6);
        pt[6 as usize] = Some(0);
        pt[10 as usize] = Some(3);
        pt[3 as usize] = Some(10);

        let parent = region(0, 10, None);
        let child = region(4, 8, Some(0)); // ri=4, inside [0,6]

        assert_eq!(
            location_status(&child, Some(&parent), &pt),
            LocationStatus::InBand
        );
    }

    #[test]
    fn test_out_band_nested_borders() {
        // parent (0,21), pt[0]=9, pt[21]=12 -> nested borders (p_bi < p_bj),
        // gap (9,12) is out-band
        let mut pt = PairTable::new(22);
        pt[0 as usize] = Some(9);
        pt[9 as usize] = Some(0);
        pt[21 as usize] = Some(12);
        pt[12 as usize] = Some(21);

        let parent = region(0, 21, None);

        let in_gap = region(10, 11, Some(0)); // ri=10, in gap (9,12)
        assert_eq!(
            location_status(&in_gap, Some(&parent), &pt),
            LocationStatus::OutBand
        );

        let left_band = region(3, 8, Some(0)); // ri=3, within [0,9]
        assert_eq!(
            location_status(&left_band, Some(&parent), &pt),
            LocationStatus::InBand
        );

        let right_band = region(15, 18, Some(0)); // ri=15, within [12,21]
        assert_eq!(
            location_status(&right_band, Some(&parent), &pt),
            LocationStatus::InBand
        );
    }
}

#[cfg(test)]
mod collect_bands_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_non_pseudo_region_has_no_bands() {
        let pt = PairTable::try_from("(((...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let region = &tree.nodes[tree.root_children[0]];
        assert!(collect_bands(&tree, region, &pt).is_empty());
    }

    #[test]
    fn test_h_type_pseudoknot_two_bands() {
        // ((([[[)))...]]] -> single region (0,14), pseudoknotted, no children.
        // Round chain (0,8)/(1,7)/(2,6) and square chain (3,14)/(4,13)/(5,12)
        // each form a 3-rung band.
        let pt = PairTable::try_from("((([[[)))...]]]").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let region = &tree.nodes[tree.root_children[0]];
        assert_eq!((region.i, region.j), (0, 14));
        assert!(region.children.is_empty());

        let bands = collect_bands(&tree, region, &pt);
        assert_eq!(bands, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }
}

#[cfg(test)]
mod nested_pairs_tests {
    use super::*;
    use std::convert::TryFrom;

    #[test]
    fn test_multiloop_children_included() {
        // (((..((...))...((...))...)))
        // -> region (2,25) is a multiloop with children (5,11) and (15,21)
        let pt = PairTable::try_from("(((..((...))...((...))...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let r0 = tree.root_children[0]; // (0,27)
        let r1 = tree.nodes[r0].children[0]; // (1,26)
        let r2 = tree.nodes[r1].children[0]; // (2,25)

        let region = &tree.nodes[r2];
        assert_eq!((region.i, region.j), (2, 25));
        assert_eq!(region.children.len(), 2);

        let pairs = nested_pairs(&tree, &region.children, &pt, 2, 25);
        assert_eq!(pairs, vec![(5, 11), (15, 21)]);
    }

    #[test]
    fn test_strict_bounds_exclude_boundary_children() {
        let pt = PairTable::try_from("(((..((...))...((...))...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let r0 = tree.root_children[0];
        let r1 = tree.nodes[r0].children[0];
        let r2 = tree.nodes[r1].children[0];
        let region = &tree.nodes[r2];

        // left == first child's i, right == second child's j -> both excluded
        let pairs = nested_pairs(&tree, &region.children, &pt, 5, 21);
        assert!(pairs.is_empty());
    }

    #[test]
    fn test_no_children() {
        let pt = PairTable::try_from("(((...)))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let r0 = tree.root_children[0];
        let r1 = tree.nodes[r0].children[0];
        let r2 = tree.nodes[r1].children[0]; // (2,6), no children

        let region = &tree.nodes[r2];
        let pairs = nested_pairs(&tree, &region.children, &pt, region.i, region.j);
        assert!(pairs.is_empty());
    }
}
