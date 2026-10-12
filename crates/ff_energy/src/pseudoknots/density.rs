//! Structure density: for each pseudoloop, the largest number of its bands
//! crossing any one position; for a structure, the maximum over its pseudoloops.
//!
//! Band count is not density — a kissing hairpin has three bands and density 2,
//! and a density-2 structure may have arbitrarily many bands.
//!
//! Purely topological: no energy model is consulted in either direction.
//!
//! Zhao, Y. S. (2005). *Efficient Algorithm for RNA Pseudoknotted Secondary
//! Structure Prediction Given a Pseudoknot Free Structure.* MSc thesis, UBC.
//! Jabbari, Condon & Zhao (2008), J Comput Biol 15:139–163, §2.3 — same
//! definition verbatim, and the source Andronescu cites.
//! Andronescu et al. (2010), RNA 16:26–42 — source of "complex" = density >= 3.
//! Their gloss "the maximum number of mutually overlapping bands", read
//! strictly as pairwise crossing (a max clique), is a different and smaller
//! quantity; it is not the definition. Rastegari & Condon (2007) name
//! `abcdcadb` (HDV, `([([))]]`) as *not* density 2, which that reading would
//! wrongly admit. Do not reintroduce it.

use ff_structure::PairTable;

use crate::pseudoknots::{RegionTree, build_closed_regions_tree, collect_bands, is_pseudo};

/// Our convention; the papers leave it unstated. Scale: 1 nested, 2 H-type or
/// kissing, >= 3 complex.
const PSEUDOKNOT_FREE_DENSITY: usize = 1;

/// Andronescu et al.'s "complex" threshold.
const COMPLEX_DENSITY: usize = 3;

/// Density of one pseudoloop, together with the band geometry it was derived
/// from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PseudoloopDensity {
    pub region: (usize, usize),
    /// Not the density under either definition.
    pub n_bands: usize,
    /// Zhao's eq. 1: most bands spanning one position.
    pub zhao: usize,
    /// R&C as printed: most bands whose *left arm* strictly contains one
    /// position. Not the definition; see the module docs.
    pub rastegari: usize,
    pub band_lengths: Vec<usize>,
}

/// `collect_bands` orders left-arm positions outer -> inner.
fn band_span(chain: &[usize], pt: &PairTable) -> (usize, usize) {
    let i = chain[0];
    let j = pt[i].expect("band rung must be paired") as usize;
    (i, j)
}

/// Largest number of bands of one crossing-graph component covering a position.
/// Zhao (2005) eq. 1: the most bands spanning any one position.
///
/// Coverage rises by one at each band start and falls at each end, so its
/// maximum is attained at some start; the other positions need not be examined.
/// O(B^2) in the band count, independent of sequence length.
fn max_coverage(spans: &[(usize, usize)]) -> usize {
    spans
        .iter()
        .map(|&(cut, _)| spans.iter().filter(|&&(i, j)| i <= cut && cut < j).count())
        .max()
        .unwrap_or(0)
}

/// `None` for a region that is not pseudoknotted.
/// R&C (2007) as printed: `i' < k < j'`, where `j'` is the band's *inner* left
/// endpoint, so k is confined to the left arm. Kept to make the discrepancy
/// reproducible, not because it is the definition.
fn max_left_arm(arms: &[(usize, usize)]) -> usize {
    // open intervals, so coverage rises at i+1, not at i
    arms.iter()
        .map(|&(i, _)| i + 1)
        .map(|k| arms.iter().filter(|&&(i, j)| i < k && k < j).count())
        .max()
        .unwrap_or(0)
}

fn region_density(tree: &RegionTree, idx: usize, pt: &PairTable) -> Option<PseudoloopDensity> {
    let region = &tree.nodes[idx];
    if !is_pseudo(region, pt) {
        return None;
    }

    let bands = collect_bands(tree, region, pt);
    if bands.is_empty() {
        return None;
    }

    let band_spans: Vec<(usize, usize)> = bands.iter().map(|c| band_span(c, pt)).collect();
    let left_arms: Vec<(usize, usize)> = bands.iter().map(|c| (c[0], c[c.len() - 1])).collect();
    let band_lengths: Vec<usize> = bands.iter().map(|c| c.len()).collect();

    Some(PseudoloopDensity {
        region: (region.i, region.j),
        n_bands: bands.len(),
        zhao: max_coverage(&band_spans),
        rastegari: max_left_arm(&left_arms),
        band_lengths,
    })
}

fn visit(tree: &RegionTree, idx: usize, pt: &PairTable, out: &mut Vec<PseudoloopDensity>) {
    if let Some(d) = region_density(tree, idx, pt) {
        out.push(d);
    }
    for &child in &tree.nodes[idx].children {
        visit(tree, child, pt, out);
    }
}

/// Density of every pseudoloop, in tree order; empty if pseudoknot-free.
///
/// Walks `root_children`, not `tree.nodes`: the arena also holds discarded
/// construction candidates that look like real top-level regions.
pub fn density_annotation(pt: &PairTable) -> Vec<PseudoloopDensity> {
    let tree = build_closed_regions_tree(pt);
    let mut out = Vec::new();
    for &idx in &tree.root_children {
        visit(&tree, idx, pt, &mut out);
    }
    out
}

/// Density, as Zhao (2005) eq. 1 defines it. This is *the* density.
pub fn zhao_density(pt: &PairTable) -> usize {
    density_annotation(pt)
        .iter()
        .map(|d| d.zhao)
        .max()
        .unwrap_or(PSEUDOKNOT_FREE_DENSITY)
}

/// What R&C's printed condition `i' < k < j'` yields. Provided so the
/// discrepancy can be reproduced; it is a transcription slip, not a rival
/// definition, and it disagrees with [`zhao_density`] on real structures.
pub fn rastegari_density(pt: &PairTable) -> usize {
    density_annotation(pt)
        .iter()
        .map(|d| d.rastegari)
        .max()
        .unwrap_or(PSEUDOKNOT_FREE_DENSITY)
}

/// "Complex" in the sense of Andronescu et al.; the complement is Zhao's
/// density-2 class.
pub fn is_complex(pt: &PairTable) -> bool {
    zhao_density(pt) >= COMPLEX_DENSITY
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pseudoknots::extended_dot_bracket_to_pair_table;
    use ff_structure::ExtendedDotBracketVec;

    fn pt_of(s: &str) -> PairTable {
        let edb = ExtendedDotBracketVec::try_from(s).unwrap();
        extended_dot_bracket_to_pair_table(&edb).unwrap()
    }

    /// The four topologies from `enumerate::band_counts_across_pseudoknot_topologies`,
    /// now asserting the density that the band count does *not* give.
    #[test]
    fn density_differs_from_band_count() {
        let cases = [
            // structure,                                       bands, density
            ("(((...[[[...)))...]]]", 2, 2),
            ("((((...[[[[...))))....((((...]]]]...))))", 3, 2),
            ("(((...[[[...{{{...)))...]]]...}}}", 3, 3),
            (
                "(((...[[[...{{{...)))...]]]...<<<...}}}...>>>",
                4,
                3,
            ),
        ];

        for (s, bands, expected) in cases {
            let pt = pt_of(s);
            let ann = density_annotation(&pt);
            let outer = ann.first().expect("one pseudoloop");
            assert_eq!(outer.n_bands, bands, "band count for {s}");
            assert_eq!(outer.zhao, expected, "density for {s}");
        }
    }

    /// Andronescu et al.: "Structures that have a density of, at most, 2
    /// include H-type pseudoknots, kissing hairpins, and pseudoknot-free
    /// structures."
    #[test]
    fn density_at_most_two_class() {
        assert_eq!(zhao_density(&pt_of("(((...[[[...)))...]]]")), 2);
        assert_eq!(
            zhao_density(&pt_of("((((...[[[[...))))....((((...]]]]...))))")),
            2
        );
        for free in ["(((...)))", "((...))..((...))", "((((...))..((...))))", "....."] {
            assert_eq!(
                zhao_density(&pt_of(free)),
                PSEUDOKNOT_FREE_DENSITY,
                "pseudoknot-free: {free}"
            );
        }
    }

    #[test]
    fn complex_is_density_three_or_more() {
        assert!(!is_complex(&pt_of("(((...[[[...)))...]]]")));
        assert!(!is_complex(&pt_of("((((...[[[[...))))....((((...]]]]...))))")));
        assert!(is_complex(&pt_of("(((...[[[...{{{...)))...]]]...}}}")));
    }

    /// Nested pseudoknots are counted per pseudoloop, not accumulated.
    #[test]
    fn nested_bands_do_not_raise_density() {
        // Two separate H-type pseudoknots, one nested inside the other.
        let pt = pt_of("((( [[[ ((( [[[ ... ))) ]]] ))) ]]]".replace(' ', "").as_str());
        let ann = density_annotation(&pt);
        assert!(ann.len() >= 2, "expected nested pseudoloops, got {ann:?}");
        for d in &ann {
            assert_eq!(d.zhao, 2, "each H-type stays at density 2: {d:?}");
        }
        assert_eq!(zhao_density(&pt), 2);
    }

    /// Counts bands spanning a position, crossing or not.
    #[test]
    fn coverage_helper() {
        assert_eq!(max_coverage(&[]), 0);
        assert_eq!(max_coverage(&[(0, 10)]), 1);
        assert_eq!(max_coverage(&[(0, 4), (6, 10)]), 1); // disjoint
        assert_eq!(max_coverage(&[(0, 10), (5, 20), (15, 30)]), 2); // chain a-b-c
        assert_eq!(max_coverage(&[(0, 10), (5, 15), (7, 20)]), 3); // mutually crossing
    }

    /// Band length supports a "needs at least one stack" filter without
    /// re-deriving the bands.
    #[test]
    fn band_lengths_distinguish_stacked_from_lone_pairs() {
        let stacked = density_annotation(&pt_of("(((...[[[...)))...]]]"));
        assert_eq!(stacked[0].band_lengths, vec![3, 3]);

        let lone = density_annotation(&pt_of("(...[...)...]"));
        assert_eq!(lone[0].band_lengths, vec![1, 1]);
    }

    /// The two readings disagree, and R&C's is usually vacuous. `((([)))]`
    /// has bands A=(0,6)(1,5)(2,4) and B=(3,7): k=3 is spanned by both, but
    /// the only k strictly inside a left arm is k=1, where only A spans.
    #[test]
    fn readings_disagree() {
        let pt = pt_of("((([)))]");
        assert_eq!(zhao_density(&pt), 2);
        assert_eq!(rastegari_density(&pt), 1);

        // single-rung bands leave the left-arm condition unsatisfiable
        let hdv = pt_of("([([))]]");
        assert_eq!(zhao_density(&hdv), 4);
        assert_eq!(rastegari_density(&hdv), 0);
    }

    /// R&C's second non-density-2 structure: pattern `abcabc`, "found at the
    /// first ribosome initiation site in the E. coli mRNA" (Gluick et al. 1997).
    /// Arcs a=(0,3) b=(1,4) c=(2,5) by their pattern convention (§5); all three
    /// cross, so density 3 under either reading.
    #[test]
    fn abcabc_is_not_density_2() {
        for s in ["([{)]}", "((([[[{{{)))]]]}}}"] {
            let pt = pt_of(s);
            let d = density_annotation(&pt).remove(0);
            assert_eq!(d.n_bands, 3, "{s}");
            assert_eq!(d.zhao, 3, "three mutually crossing bands: {s}");
            assert!(is_complex(&pt), "R&C: not a density 2 structure: {s}");
        }
    }

    /// R&C name pattern `abcdcadb` (HDV) as "not a density 2 structure". Arcs
    /// a=(0,5) b=(1,7) c=(2,4) d=(3,6) by their pattern convention (§5); ab, ad
    /// and cd form a tree, so a max-clique reading gives 2 and wrongly admits
    /// it. Counting bands over a position gives 4, matching the paper.
    #[test]
    fn hdv_is_not_density_2() {
        let pt = pt_of("([([))]]");

        let d = density_annotation(&pt).remove(0);
        assert_eq!(d.n_bands, 4);
        assert_eq!(d.zhao, 4, "all four bands span position 3.5");

        assert_eq!(zhao_density(&pt), 4);
        assert!(is_complex(&pt), "R&C: not a density 2 structure");

        // Two bracket families, yet outside the density-2 class.
        assert!(!"([([))]]".contains(['{', '<']));
    }
}


