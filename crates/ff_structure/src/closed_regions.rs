//! Closed regions and bands of (pseudoknotted) secondary structures.
//!
//! A closed region is an interval [i, j] such that every base pair with one
//! end in it has both ends in it, and no smaller interval starting at i has
//! this property. It is nested if i pairs with j, and pseudoknotted
//! otherwise. The pairs of a pseudoknotted region that are not part of a
//! region inside it form its bands: chains of pairs (outer -> inner) where
//! consecutive pairs enclose nothing but other regions.
//!
//! Strand breaks: two consecutive pairs of a band are only chained if no
//! strand break lies between them (outside of skipped regions). A strand
//! break there means the two pairs are not connected by a loop, so they form
//! two bands. A strand wrap counts as a strand break here: in the order
//! of positions it separates strands.
//!
//! These bands follow Rastegari & Condon: a band runs across nested regions
//! between its pairs. The strained stems of `ff_energy::Pseudoknot` do not:
//! a stem ends at any loop other than a stack, bulge or interior loop.
//! Used for the pseudoknot density (`ff-eval -v`).
//!
//! TODO: compare with the literature: whether a band chains across a
//! pseudoknotted region inside it, and the density (Zhao).

use crate::MultiStruct;
use crate::MultiPairTable;

/// A closed region, see [`ClosedRegions`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClosedRegion {
    /// The first and the last index of the region.
    pub span: (usize, usize),
    /// True if the first and the last base of the region pair with each other.
    pub nested: bool,
}

/// The closed regions of a structure, built in one pass: a stack of open
/// regions, merged when a pair crosses them. Borrows the table it was built
/// from.
///
/// A region is a closed region: an interval [i, j] such that every pair
/// with one end in it has both ends in it, and no smaller interval starting
/// at i has this property. It is nested if i pairs with j, pseudoknotted
/// otherwise. Every region starts at an opening base, but not every opening
/// base starts a region: inside a pseudoknotted region, only the first base
/// does. Regions are found with [`at`](Self::at) by their first index.
///
/// For a circular strand, regions and bands are those of the drawing, cut at
/// the strand wrap: a helix through the cut is two helices there, so they
/// change when the circle is rotated. Whether a structure has a pseudoknot
/// does not depend on the cut.
///
/// A structure analysis (closed regions, bands, density); the energy
/// evaluation of `ff_energy` does not use it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ClosedRegions<'a> {
    table: &'a MultiPairTable,
    /// The region starting at each index, if any.
    regions: Vec<Option<ClosedRegion>>,
}

impl<'a> From<&'a MultiPairTable> for ClosedRegions<'a> {
    fn from(table: &'a MultiPairTable) -> Self {
        let n = table.len();
        let mut regions = vec![None; n];
        // Open regions as (first, last) index, last extended while merging.
        let mut stack: Vec<(usize, usize)> = Vec::new();

        for k in 0..n {
            if let MultiStruct::Paired(p) = table.get(k) {
                if k < p {
                    stack.push((k, p));
                } else {
                    // Pair (p, k) closes: every open region that starts
                    // after p crosses it, merge them into the region below.
                    let mut last = k;
                    while let Some(&(first, end)) = stack.last() && first > p {
                        last = last.max(end);
                        stack.pop();
                    }
                    if let Some(top) = stack.last_mut() {
                        top.1 = top.1.max(last);
                    }
                }
            }
            if let Some(&(first, last)) = stack.last() && last == k {
                stack.pop();
                let nested = table.get(first) == MultiStruct::Paired(last);
                regions[first] = Some(ClosedRegion {
                    span: (first, last),
                    nested,
                });
            }
        }
        debug_assert!(stack.is_empty(), "Unbalanced pair table.");
        Self { table, regions }
    }
}

impl ClosedRegions<'_> {

    /// The region starting at index p, if any.
    pub fn at(&self, p: usize) -> Option<ClosedRegion> {
        self.regions[p]
    }

    /// All closed regions, sorted by their first index.
    pub fn iter(&self) -> impl Iterator<Item = ClosedRegion> + '_ {
        self.regions.iter().flatten().copied()
    }

    /// The bands of a pseudoknotted region, each as its pairs outer -> inner,
    /// sorted by the first index of the outer pair. All regions inside 
    /// `region` are skipped, nested or pseudoknotted. Used for the density.
    fn region_bands(&self, region: ClosedRegion) -> Vec<Vec<(usize, usize)>> {
        debug_assert!(!region.nested, "Bands of a nested region.");
        self.table.band_walk(region.span, |k| self.regions[k]
            .filter(|_| k != region.span.0)
            .map(|child| child.span.1))
    }

    /// The density of the structure: the largest number of bands of one
    /// pseudoknotted region that span a common position. 1 for nested
    /// structures, 2 for H-types and kissing hairpins. This is not the largest
    /// number of pairwise crossing bands.
    ///
    /// Bands are split at strand breaks, so in multi-stranded structures the
    /// density depends on the order of strands.
    pub fn density(&self) -> usize {
        self.iter()
            .filter(|r| !r.nested)
            .map(|r| {
                let spans: Vec<(usize, usize)> = self.region_bands(r).iter()
                    .map(|band| band[0])
                    .collect();
                // Coverage only increases at the start of a band.
                spans.iter()
                    .map(|&(k, _)| spans.iter().filter(|&&(p, q)| p <= k && k < q).count())
                    .max()
                    .unwrap_or(0)
            })
            .max()
            .unwrap_or(1)
    }
}

impl MultiPairTable {
    /// The bands in `span`, each as its pairs outer -> inner, sorted by the
    /// first index of the outer pair. `child_end(k)` returns the last index
    /// of a region starting at k that is skipped, if any. Consecutive pairs
    /// are not chained across a strand break outside of skipped regions.
    ///
    /// Collects the paired bases outside of skipped regions, then chains
    /// each pair with the next paired bases inward on both sides, as long as
    /// they pair with each other.
    fn band_walk(
        &self,
        (i, j): (usize, usize),
        child_end: impl Fn(usize) -> Option<usize>,
    ) -> Vec<Vec<(usize, usize)>> {
        let partner = |k: usize| match self.get(k) {
            MultiStruct::Paired(p) => p,
            _ => unreachable!("Band bases are paired."),
        };

        // The paired bases of the region outside of nested regions (R&C's BL), 
        // and the strand breaks among them, both 5'->3'.
        let mut paired: Vec<usize> = Vec::new();
        let mut breaks: Vec<usize> = Vec::new();
        let mut k = i;
        while k <= j {
            match child_end(k) {
                Some(last) => k = last,
                None => match self.get(k) {
                    MultiStruct::Paired(_) => paired.push(k),
                    MultiStruct::StrandBreak | MultiStruct::StrandWrap => breaks.push(k),
                    MultiStruct::Unpaired => (),
                },
            }
            k += 1;
        }
        let no_break = |a: usize, b: usize| !breaks.iter().any(|&c| a < c && c < b);

        let mut bands = Vec::new();
        let mut r = 0;
        while r < paired.len() {
            let b = partner(paired[r]);
            if paired[r] > b { // the 3' base of a pair
                r += 1;
                continue;
            }
            let (mut ra, mut rb) = (r, paired.binary_search(&b).expect("A band base."));
            let mut band = vec![(paired[ra], b)];
            // Chain the next pair inward while it directly follows on both sides.
            while ra + 1 < rb {
                let (c, d) = (paired[ra + 1], paired[rb - 1]);
                if partner(c) != d || !no_break(paired[ra], c) || !no_break(d, paired[rb]) {
                    break;
                }
                band.push((c, d));
                (ra, rb) = (ra + 1, rb - 1);
            }
            bands.push(band);
            r = ra + 1;
        }
        bands
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(s: &str) -> MultiPairTable {
        MultiPairTable::try_from(s).unwrap()
    }

    fn spans(r: &ClosedRegions) -> Vec<(usize, usize, bool)> {
        r.iter().map(|r| (r.span.0, r.span.1, r.nested)).collect()
    }

    #[test]
    fn nested_regions() {
        let t = table("(((...)))");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![(0, 8, true), (1, 7, true), (2, 6, true)]);
        assert!(t.is_nested());

        let t = table("((..))((..))");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![(0, 5, true), (1, 4, true), (6, 11, true), (7, 10, true)]);

        let t = table(".....");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![]);
    }

    #[test]
    fn h_type_is_one_region() {
        let t = table("((([[[)))...]]]");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![(0, 14, false)]);
        assert_eq!(r.at(1), None); // a band pair, not a region
        assert!(!t.is_nested());
    }

    #[test]
    fn turnip_yellow_mosaic_virus() {
        let t = table("(((((.......)))))(((....[[[[[)))...]]]]]....");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![
            (0, 16, true), (1, 15, true), (2, 14, true), (3, 13, true), (4, 12, true),
            (17, 39, false),
        ]);
    }

    #[test]
    fn nested_regions_inside_a_pseudoknot() {
        // (3,7) lies in a band region, (27,31) between band regions.
        let t = table("((.(...)..((..[[[..))..))..(...)...]]]");
        let r = ClosedRegions::from(&t);
        assert_eq!(spans(&r), vec![(0, 37, false), (3, 7, true), (27, 31, true)]);
    }

    #[test]
    fn bands_h_type() {
        let t = table("((([[[)))...]]]");
        let r = ClosedRegions::from(&t);
        let pk = r.at(0).unwrap();
        assert_eq!(r.region_bands(pk), vec![
            vec![(0, 8), (1, 7), (2, 6)],
            vec![(3, 14), (4, 13), (5, 12)],
        ]);
    }

    #[test]
    fn bands_skip_nested_regions() {
        // Band A continues across the multiloop with child (3,7).
        let t = table("((.(...)..((..[[[..))..))..(...)...]]]");
        let r = ClosedRegions::from(&t);
        let pk = r.at(0).unwrap();
        assert_eq!(r.region_bands(pk), vec![
            vec![(0, 24), (1, 23), (10, 20), (11, 19)],
            vec![(14, 37), (15, 36), (16, 35)],
        ]);
    }

    #[test]
    fn bands_split_at_strand_breaks() {
        // [[[[...((((+((((...]]]]+))))...))))  (C+A+B): Z and X would form
        // one band, but a strand break lies between them.
        let t = table("[[[[...((((+((((...]]]]+))))...))))");
        let r = ClosedRegions::from(&t);
        let pk = r.at(0).unwrap();
        let outer = |b: &Vec<(usize, usize)>| b[0];
        let bands = r.region_bands(pk);
        assert_eq!(bands.iter().map(outer).collect::<Vec<_>>(), vec![(0, 22), (7, 34), (12, 27)]);
        assert!(bands.iter().all(|b| b.len() == 4));

        // Without the strand break, they form one band.
        let t = table("[[[[...((((.((((...]]]]+))))...))))");
        let r = ClosedRegions::from(&t);
        let bands = r.region_bands(r.at(0).unwrap());
        assert_eq!(bands.iter().map(outer).collect::<Vec<_>>(), vec![(0, 22), (7, 34)]);
        assert_eq!(bands[1].len(), 8);
    }

    #[test]
    fn hidden_strand_breaks_do_not_split() {
        // The break lies inside the nested region (11,13) between band pairs.
        let t = table("[[[[...((((.(+)((((...]]]].))))...))))");
        let r = ClosedRegions::from(&t);
        let bands = r.region_bands(r.at(0).unwrap());
        assert_eq!(bands.len(), 2);
        assert_eq!(bands[1].len(), 8);
    }

    #[test]
    fn density() {
        for (s, n_bands, density) in [
            ("(((...[[[...)))...]]]", 2, 2),
            ("((((...[[[[...))))....((((...]]]]...))))", 3, 2), // kissing hairpin
            ("(((...[[[...{{{...)))...]]]...}}}", 3, 3),
            ("(((...[[[...{{{...)))...]]]...<<<...}}}...>>>", 4, 3),
            ("([{)]}", 3, 3),
            ("((([[[{{{)))]]]}}}", 3, 3),
            ("([([))]]", 4, 4), // HDV-like: not density 2
            ("((([)))]", 2, 2),
        ] {
            let t = table(s);
        let r = ClosedRegions::from(&t);
            let pk = r.iter().find(|r| !r.nested).unwrap();
            assert_eq!(r.region_bands(pk).len(), n_bands, "bands of {s}");
            assert_eq!(r.density(), density, "density of {s}");
        }
        for s in ["(((...)))", "((...))..((...))", "((((...))..((...))))", "....."] {
            let t = table(s);
        let r = ClosedRegions::from(&t);
            assert_eq!(r.density(), 1, "density of {s}");
        }
    }

    #[test]
    fn knots_in_stretches_are_regions() {
        // A knot in a stretch of another: two regions, each with its bands.
        let t = table("(((<<<[[[{{{...]]]}}})))>>>");
        let r = ClosedRegions::from(&t);
        let outer = |bands: Vec<Vec<(usize, usize)>>| bands.iter().map(|b| b[0]).collect::<Vec<_>>();
        assert_eq!(outer(r.region_bands(r.at(0).unwrap())), vec![(0, 23), (3, 26)]);
    }

    #[test]
    fn nested_pseudoknots_do_not_add_up() {
        // Two H-types, one inside the other: density 2.
        let t = table("(((<<<[[[{{{...]]]}}})))>>>");
        let r = ClosedRegions::from(&t);
        assert_eq!(r.iter().filter(|r| !r.nested).count(), 2);
        assert_eq!(r.density(), 2);
    }
}
