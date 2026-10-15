//! Energy parameters for pseudoknotted structures.
//!
//! All energy values are in **dcal/mol** (1/100 kcal/mol), consistent with
//! the rest of `ff_energy`.
//!
//! The full Dirks-Pierce model has 11 features (Andronescu et al. 2010, Table 6).
//!
//! The Pseudoloop energy formula:
//!
//! ```text
//! E = init(ctx) + pb × n_bands + pup × n_pup + pps × n_nested
//! ```
//!
//! SpanBand multiloop formula:
//!
//! ```text
//! E = ap + bp × n_branches + cp × n_unpaired
//! ```
//!
//! Stacks and interior loops **that span a band** (SpanBand loops in fuzzyfold) use
//! scaled Turner energies rather than the raw Turner value:
//! - stacked pair: `round(e_stp × turner_stack)`
//! - interior loop / bulge: `round(e_intp × turner_interior)`
//!
//! This is the multiplicative penalty from Andronescu et al. (2010) §Materials and Methods.
//! It makes the full DP energy function quadratic in the parameters (the SpanBand count
//! is multiplied by both `e_stp`/`e_intp` and the underlying Turner parameter), but each
//! loop is evaluated independently so the implementation remains linear in structure size.

/// Full Dirks-Pierce (DP) energy parameter set for pseudoknots.
///
/// Covers all 11 features described in Andronescu, Pop & Condon (2010), Table 6.
/// Two static constants are provided:
/// - [`RNA_DP03`]: the original Dirks & Pierce 2003 values.
/// - [`RNA_DP09`]: Andronescu 2010 retrained values.
///
/// The DP params are independent of the nearest-neighbor model, so any combination
/// evaluates. Only one *trained* pairing is reachable from this crate:
/// DP09 = MT09 + dp09 (Andronescu et al. 2010, Table 3 row 9). The paper's other
/// named set, DP03 = MT99 + dp03 (row 2), pairs dp03 with Mathews et al. (1999),
/// which this crate does not carry — so DP03 cannot be expressed here. Turner 2004
/// is a different nearest-neighbor model, not a stand-in for MT99.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DPParams {
    /// Initiation when the pseudoloop is in an exterior context (dcal/mol).
    pub init_external: i32,
    /// Initiation when the pseudoloop is nested in a standard multiloop (dcal/mol).
    pub init_multiloop: i32,
    /// Initiation when the pseudoloop is nested in another pseudoloop (dcal/mol).
    pub init_pseudoloop: i32,
    /// Per-band (per crossing helix) penalty (dcal/mol).
    pub pb: i32,
    /// Per-free-linker-nucleotide penalty (dcal/mol).
    /// Applied to `n_pup`: the total count of unpaired bases across all 2n−1
    /// inter-arm gaps that belong to the pseudoloop itself (i.e. not inside
    /// any nested child region).
    pub pup: i32,
    /// Per-nested-closed-region inside the pseudoloop (dcal/mol).
    pub pps: i32,
    /// SpanBand multiloop: initiation penalty (dcal/mol).
    pub ap: i32,
    /// SpanBand multiloop: per-helix penalty (dcal/mol).
    pub bp: i32,
    /// SpanBand multiloop: per-unpaired-base penalty (dcal/mol).
    pub cp: i32,
    /// Multiplicative scale factor for Turner stacking energies that span a band.
    /// Applied as `round(e_stp × turner_stack)` for each SpanBand Stack loop.
    pub e_stp: f64,
    /// Multiplicative scale factor for Turner interior-loop / bulge energies that span a band.
    /// Applied as `round(e_intp × turner_interior)` for each SpanBand Interior/Bulge loop.
    pub e_intp: f64,
}

/// Initial Dirks-Pierce 2003 pseudoknot parameters for RNA at 37 °C.
///
/// Source: Andronescu, Pop & Condon, *RNA* **16** (2010), Table 6, column dp03.
/// These are the parameter values used by the original HotKnots / NUPACK software.
/// Published as DP03 = MT99 + dp03 (Table 3 row 2). MT99 is Mathews et al.
/// (1999); no MT99 table exists in this crate, so DP03 itself is out of reach.
pub static RNA_DP03: DPParams = DPParams {
    init_external: 960,    //  9.60 kcal/mol
    init_multiloop: 1500,  // 15.00 kcal/mol
    init_pseudoloop: 1500, // 15.00 kcal/mol
    pb: 20,                //  0.20 kcal/mol per band
    pup: 10,               //  0.10 kcal/mol per unpaired base
    pps: 10,               //  0.10 kcal/mol per nested closed region
    ap: 340,               //  3.40 kcal/mol SpanBand multiloop initiation
    bp: 40,                //  0.40 kcal/mol SpanBand per helix
    cp: 0,                 //  0.00 kcal/mol SpanBand per unpaired base
    e_stp: 0.83,
    e_intp: 0.83,
};

/// Trained Dirks-Pierce 2009 pseudoknot parameters for RNA at 37 °C.
///
/// Source: Andronescu, Pop & Condon, *RNA* **16** (2010), Table 6, column dp09.
/// These are the best-performing DP parameters (F-measure 79% vs. 68% for dp03).
/// They are identical to the HotKnots v2 and Spark `PK_globals.hh` values.
/// Published as DP09 = MT09 + dp09 ([`RNA_MT09`], Table 3 row 9) — the one
/// trained pairing this crate can reproduce.
///
/// Key differences from dp03:
/// - External initiation flips from +960 → −138 (exterior pseudoknots become rewarded)
/// - Per-band penalty rises from 20 → 246 (discourages over-complex pseudoknots)
///
/// NOTE on `init_external`. Table 6 of the published paper prints this value as
/// `1.38`, without a sign. The authoritative source is negative:
///
/// ```text
/// // HotKnots v2 / Spark, src/PK_globals.hh
/// int PS_penalty = -138;  // 960;  // exterior pseudoloop initiation penalty
/// ```
///
/// The `// 960` trailing comment is the dp03 value it replaced, which matches
/// Table 6's dp03 column (9.60) — so the table and the source agree everywhere
/// except for this dropped minus sign. All eleven values here were verified
/// field-by-field against `PK_globals.hh`.
pub static RNA_DP09: DPParams = DPParams {
    init_external: -138,  // -1.38 kcal/mol (slight reward — exterior PK is favorable)
    init_multiloop: 1007, // 10.07 kcal/mol
    init_pseudoloop: 1500, // 15.00 kcal/mol
    pb: 246,              //  2.46 kcal/mol per band
    pup: 6,               //  0.06 kcal/mol per unpaired base
    pps: 96,              //  0.96 kcal/mol per nested closed region
    ap: 341,              //  3.41 kcal/mol SpanBand multiloop initiation
    bp: 56,               //  0.56 kcal/mol SpanBand per helix
    cp: 12,               //  0.12 kcal/mol SpanBand per unpaired base
    e_stp: 0.89,
    e_intp: 0.74,
};

// Cao & Chen (CC09) is not implemented: its energy function is non-linear
// (min(coaxial, dangle)) and its loop entropies are ΔS°, neither of which fits
// the linear loop-sum here. It also does not outperform DP09.
