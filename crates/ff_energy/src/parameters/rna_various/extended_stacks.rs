//! # Pseudouridine Parameters
//!
//! **Authors:** Graham A. Hudson, Richard J. Bloomingdale, and Brent M. Znosko  
//! **Title:** *Thermodynamic contribution and nearest-neighbor parameters of
//! pseudouridine–adenosine base pairs in oligoribonucleotides*  
//! **Journal:** RNA 19:1474–1482  
//! **Year:** 2013  
//! **DOI:** 10.1261/rna.039610.113
//!
//! ## Description
//!
//! This module provides stacking parameters including **P** for
//! pseudouridine (Ψ).
//!
//! In the implementation:
//!
//! - All unspecified pseudouridine interactions are treated as **U**.
//! - NOTE: The corrected terminal **PU** end penalty is NOT applied.
//!
//! ```rust
//! let terminal_pu_en37  =  31;
//! let terminal_pu_enth  = -204;
//! ```
//!
//! ## Implementation Notes
//!
//! The handling is slightly subtle because the terminal mismatch
//! contribution is already baked into selected nearest-neighbor
//! parameters. Care must be taken to avoid double-counting or
//! inconsistent corrections when applying PU-specific penalties.
//!
use crate::parameters::ExtendedStackParams;

pub static STACK_EN37: ExtendedStackParams = [
    /* [cl] [i]:   AU     UA     CG     GC     GU     UG     AP     PA */
    /* [AU] */ [ -110,   -90,  -210,  -220,  -140,   -60,  -274,  -280],
    /* [UA] */ [  -90,  -130,  -210,  -240,  -130,  -100,  -210,  -162],
    /* [CG] */ [ -210,  -210,  -240,  -330,  -210,  -140,  -220,  -277],
    /* [GC] */ [ -220,  -240,  -330,  -340,  -250,  -150,  -249,  -329],
    /* [GU] */ [ -140,  -130,  -210,  -250,   130,   -50,  -140,  -130],
    /* [UG] */ [  -60,  -100,  -140,  -150,   -50,    30,   -60,  -100],
    /* [AP] */ [ -162,  -280,  -329,  -277,  -140,   -60,  -110,   -90],
    /* [PA] */ [ -210,  -274,  -249,  -220,  -130,  -100,   -90,  -130],
];

pub static STACK_ENTH: ExtendedStackParams = [
    /* [cl] [i]:  AU     UA     CG     GC     GU     UG     AP     PA */
    /* [AU] */ [ -940,  -680, -1050, -1140,  -880,  -320, -2694, -2208],
    /* [UA] */ [ -680,  -770, -1040, -1240, -1280,  -700, -1247, -2081],
    /* [CG] */ [-1050, -1040, -1060, -1340, -1210,  -560, -1119, -1623],
    /* [GC] */ [-1140, -1240, -1340, -1490, -1260,  -830, -1729, -2407],
    /* [GU] */ [ -880, -1280, -1210, -1260, -1460, -1350,  -880, -1280],
    /* [UG] */ [ -320,  -700,  -560,  -830, -1350,  -930,  -320,  -700],
    /* [AP] */ [-2081, -2208, -2407, -1623,  -880,  -320,  -940,  -680],
    /* [PA] */ [-1247, -2694, -1729, -1119, -1280,  -700,  -680,  -770],
];

