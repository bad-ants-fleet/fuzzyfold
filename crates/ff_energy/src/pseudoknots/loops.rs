//! loop type definitions
use std::fmt;

use crate::pseudoknots::{ClosingDescriptor, LocationStatus};

/// The seven loop categories from Rastegari & Condon's classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopType {
    Stack,
    Hairpin,
    Interior,
    Bulge,
    Multiloop,
    External,
    Pseudoloop,
}

/// Context in which a [`LoopType::Pseudoloop`] appears.
///
/// Used by the D&P energy model (features 1–3 in Andronescu et al. 2010,
/// Table 6) to select the appropriate initiation penalty:
/// - Exterior:  `init_external`   (9.60 kcal/mol dp03, 1.38 dp09)
/// - Multiloop: `init_multiloop`  (15.00 kcal/mol dp03, 10.07 dp09)
/// - Pseudoloop:`init_pseudoloop` (15.00 kcal/mol, same for both)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PseudoloopContext {
    /// Pseudoloop is at the top level (exterior).
    External,
    /// Pseudoloop is nested inside a standard multiloop.
    Multiloop,
    /// Pseudoloop is nested inside another pseudoloop.
    Pseudoloop,
}

/// Classifies the loop between a closing pair and a strictly nested inner
/// pair as a stack, bulge, or interior loop, and returns the number of
/// unpaired bases on the 5' and 3' sides.
pub fn interior_loop_type(
    closing: (usize, usize),
    inner: (usize, usize),
) -> (LoopType, usize, usize) {
    let (ci, cj) = closing;
    let (ii, ij) = inner;

    debug_assert!(
        ci < ii && ij < cj,
        "inner pair {inner:?} must be strictly nested inside closing pair {closing:?}"
    );

    let n5 = ii - ci - 1; // number of unpaired bases on 5' end
    let n3 = cj - ij - 1; // number of unpaired bases on 3' end

    let loop_type = match (n5, n3) {
        (0, 0) => LoopType::Stack,
        (0, _) | (_, 0) => LoopType::Bulge,
        _ => LoopType::Interior,
    };

    (loop_type, n5, n3)
}

/// Mirrors the Python `Loop` dataclass.
///
/// `closing` and elements of `children` use [`ClosingDescriptor`] to represent
/// either a single pair or two crossing pairs.
///
/// Fields specific to `Pseudoloop` loops:
/// - `n_pup`: total free linker nucleotides across all inter-arm gaps (the `pup`
///   term in D&P).  Computed as the sum of the 2n−1 raw positional gap spans minus
///   the span of any nested child region that falls inside each gap — matching
///   HotKnots' `pseudoEnergyDP` exactly.
/// - `n_bands`: number of crossing helices (bands) in this pseudoloop.
/// - `n_nested`: number of closed regions nested *inside* the pseudoloop
///   (distinct from the band tips listed in `children`).
/// - `pk_context`: the loop context used for the D&P initiation penalty.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Loop {
    pub loop_type: LoopType,
    pub location: LocationStatus,
    pub closing: Option<ClosingDescriptor>,
    pub inner: Option<(usize, usize)>,
    pub children: Vec<ClosingDescriptor>,
    pub unpaired_5p: usize,
    pub unpaired_3p: usize,
    /// Free linker nucleotides summed across all 2n−1 inter-arm gaps. Set on `Pseudoloop` loops.
    pub n_pup: usize,
    /// Number of bands (crossing helices). Set on `Pseudoloop` loops.
    pub n_bands: usize,
    /// Number of actual nested closed regions inside the pseudoloop.
    pub n_nested: usize,
    /// Context for context-dependent initiation penalty. Set on `Pseudoloop` loops.
    pub pk_context: Option<PseudoloopContext>,
}

impl Loop {
    pub fn new(loop_type: LoopType, location: LocationStatus) -> Self {
        Loop {
            loop_type,
            location,
            closing: None,
            inner: None,
            children: Vec::new(),
            unpaired_5p: 0,
            unpaired_3p: 0,
            n_pup: 0,
            n_bands: 0,
            n_nested: 0,
            pk_context: None,
        }
    }

    pub fn with_closing(mut self, closing: ClosingDescriptor) -> Self {
        self.closing = Some(closing);
        self
    }

    pub fn with_inner(mut self, inner: (usize, usize)) -> Self {
        self.inner = Some(inner);
        self
    }

    pub fn with_children(mut self, children: Vec<ClosingDescriptor>) -> Self {
        self.children = children;
        self
    }

    pub fn with_unpaired(mut self, n5: usize, n3: usize) -> Self {
        self.unpaired_5p = n5;
        self.unpaired_3p = n3;
        self
    }

    /// Set the total free linker nucleotide count for a pseudoloop.
    pub fn with_pup(mut self, n: usize) -> Self {
        self.n_pup = n;
        self
    }

    /// Set the number of bands (crossing helices). For `Pseudoloop` loops.
    pub fn with_bands(mut self, n: usize) -> Self {
        self.n_bands = n;
        self
    }

    /// Set the number of nested closed regions inside the pseudoloop.
    pub fn with_nested(mut self, n: usize) -> Self {
        self.n_nested = n;
        self
    }

    /// Set the context for the D&P initiation penalty. For `Pseudoloop` loops.
    pub fn with_pk_context(mut self, ctx: PseudoloopContext) -> Self {
        self.pk_context = Some(ctx);
        self
    }

    /// The closing pair(s) of this loop, flattened: 0, 1, or 2 pairs.
    /// Mirrors the Python `closing_pairs` property.
    pub fn closing_pairs(&self) -> Vec<(usize, usize)> {
        match self.closing {
            None => Vec::new(),
            Some(ClosingDescriptor::Single(p)) => vec![p],
            Some(ClosingDescriptor::Double(p1, p2)) => vec![p1, p2],
        }
    }
}

impl fmt::Display for Loop {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts: Vec<String> = vec![format!("{:?}", self.loop_type)];

        if self.location != LocationStatus::Standard {
            parts.push(format!("{:?}", self.location));
        }
        if let Some(closing) = &self.closing {
            parts.push(format!("closing={closing}"));
        }
        if let Some((i, j)) = self.inner {
            parts.push(format!("inner=({i}, {j})"));
        }
        if !self.children.is_empty() {
            let children_str = self
                .children
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            parts.push(format!("children=[{children_str}]"));
        }
        if self.unpaired_5p != 0 || self.unpaired_3p != 0 {
            parts.push(format!(
                "unpaired=({}, {})",
                self.unpaired_5p, self.unpaired_3p
            ));
        }
        if self.loop_type == LoopType::Pseudoloop {
            parts.push(format!("pup={}", self.n_pup));
            parts.push(format!("bands={}", self.n_bands));
            parts.push(format!("nested={}", self.n_nested));
            if let Some(ctx) = self.pk_context {
                parts.push(format!("ctx={ctx:?}"));
            }
        }

        write!(f, "Loop({})", parts.join(", "))
    }
}

#[cfg(test)]
mod interior_loop_tests {
    use super::*;
    use crate::{build_closed_regions_tree, closing_pairs};
    use ff_structure::PairTable;
    use std::convert::TryFrom;

    #[test]
    fn test_stack() {
        let (lt, n5, n3) = interior_loop_type((0, 8), (1, 7));
        assert_eq!(lt, LoopType::Stack);
        assert_eq!((n5, n3), (0, 0));
    }

    #[test]
    fn test_bulge_5p() {
        let (lt, n5, n3) = interior_loop_type((0, 9), (1, 7));
        assert_eq!(lt, LoopType::Bulge);
        assert_eq!((n5, n3), (0, 1));
    }

    #[test]
    fn test_bulge_3p() {
        let (lt, n5, n3) = interior_loop_type((0, 8), (2, 7));
        assert_eq!(lt, LoopType::Bulge);
        assert_eq!((n5, n3), (1, 0));
    }

    #[test]
    fn test_interior() {
        let (lt, n5, n3) = interior_loop_type((1, 12), (4, 9));
        assert_eq!(lt, LoopType::Interior);
        assert_eq!((n5, n3), (2, 2));
    }

    /// Wires interior_loop_type up to a real structure via closing_pairs.
    #[test]
    fn test_interior_loop_from_structure() {
        // ((..((..))..))
        let pt = PairTable::try_from("((..((..))..))").unwrap();
        let tree = build_closed_regions_tree(&pt);

        let r0 = tree.root_children[0]; // (0,13)
        let r1 = tree.nodes[r0].children[0]; // (1,12)
        let r2 = tree.nodes[r1].children[0]; // (4,9)

        let closing = closing_pairs(&tree.nodes[r1], &pt)[0]; // (1,12)
        let inner = closing_pairs(&tree.nodes[r2], &pt)[0]; // (4,9)

        let (lt, n5, n3) = interior_loop_type(closing, inner);
        assert_eq!(lt, LoopType::Interior);
        assert_eq!((n5, n3), (2, 2));
    }
}
