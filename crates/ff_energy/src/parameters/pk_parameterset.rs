//! Energy parameters for pseudoknotted structures.
//!
//! All energy values are in dcal/mol (1/100 kcal/mol), consistent with
//! the rest of `ff_energy`. They are given at 37 C and are not rescaled
//! with temperature.
//!
//! Which parts of a pseudoknot the terms apply to is decided by the cycles of
//! the loop graph, independent of the order of the strands (see
//! `Pseudoknot`): a stem is strained if its pairs lie on a cycle, a segment
//! of a loop if two visits of the same block enclose it.
//!
//! ```text
//! E = init     * knots (blocks that share strained segments)
//!   + stem     * strained stems
//!   + unpaired * unpaired bases in the pseudoloop
//!   + branch   * branches of the pseudoloop (nested pairs, relaxed helix ends)
//!   + sum of round(stack_scale * stack)        over stacks in strained stems
//!   + sum of round(interior_scale * interior)  over interior loops and bulges
//!                                              in strained stems
//!   + helix ends: terminal penalty, mismatch or dangles (see below)
//! ```
//!
//! Everything else is evaluated with the nearest-neighbor model, as loops
//! outside the pseudoknot: stacks and interior loops of relaxed stems, and
//! the loops on no cycle. A segment of a pseudoknot that is not strained
//! leads to a strand end: its helix ends are evaluated as in the exterior
//! loop. One duplex initiation per strand joined, for the complex.
//!
//! Every helix end of a loop of the pseudoknot gets its terminal 
//! penalty, and a mismatch if both of its flanks are in the same context 
//! (`mismatch_multi` in the pseudoloop, `mismatch_exterior` in an 
//! exterior-type stretch), else the dangle on each side (`-d2` neighbors). 
//! This includes the ends of strained stems and the children in the 
//! pseudoloop.
//!
//! TODO: compare with the Dirks & Pierce model as used in HotKnots. Its
//! counterparts of the terms are noted below; it also has an initiation in a
//! multiloop or in another pseudoloop, and parameters for a multiloop inside
//! a band, which have no counterpart here. Also compare: terminal penalties
//! and dangles of helix ends in the pseudoloop, and whether a stem should 
//! continue through large interior loops (here it does, scaled by 
//! `interior_scale`).

/// Energy parameters for pseudoknots.
///
/// The default values come from a preliminary training on the HotKnots
/// data sets (Andronescu, Pop & Condon 2010), with the nearest neighbor
/// parameters of `RNA_EXTENDED`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PKParams {
    /// Initiation, per knot: blocks of the loop graph that share strained
    /// segments.
    /// HotKnots: Ps, the initiation in an exterior context.
    pub init: i32,
    /// Per strained stem (its pairs lie on a cycle of the loop graph).
    /// HotKnots: Pb, per band.
    pub stem: i32,
    /// Per branch in a strained segment: a pair that is not on a cycle (a
    /// nested pair, a relaxed stem, an enclosing pair).
    /// HotKnots: Pps, per nested pair.
    pub branch: i32,
    /// Per unpaired base in a strained segment. HotKnots: Pup.
    pub unpaired: i32,
    /// Scale factor for stacking energies in a strained stem, applied
    /// as `round(stack_scale * stack)` to each stack. HotKnots: stP.
    pub stack_scale: f64,
    /// Scale factor for interior loop and bulge energies in a strained
    /// stem, applied as `round(interior_scale * interior)` to each loop. 
    /// HotKnots: intP.
    pub interior_scale: f64,
}

impl Default for PKParams {
    fn default() -> Self {
        PKParams {
            init: 179,             //  1.79 kcal/mol
            stem: 519,             //  5.19 kcal/mol
            branch: 164,           //  1.64 kcal/mol
            unpaired: 10,          //  0.10 kcal/mol
            stack_scale: 0.90,
            interior_scale: 1.12,
        }
    }
}
