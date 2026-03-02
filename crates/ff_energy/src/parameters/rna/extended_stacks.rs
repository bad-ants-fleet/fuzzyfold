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
//! ## Implementation Notes
//!
//! - All unspecified pseudouridine interactions are treated as **U**.
//! - [AP][AP] = [AP][AU]
//! - [AP][GU] = [AU][GU]
//! - [GP] = [GU] (not even shown in table.
//!
//!
use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::E;

pub static STACKPARAMS_EN37: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA */
    /* [AU] */ [ -110,   -90,  -210,  -220,  -140,   -60,  -280,  -274],
    /* [UA] */ [  -90,  -130,  -210,  -240,  -130,  -100,  -162,  -210],
    /* [CG] */ [ -210,  -210,  -240,  -330,  -210,  -140,  -277,  -220],
    /* [GC] */ [ -220,  -240,  -330,  -340,  -250,  -150,  -329,  -249],
    /* [GU] */ [ -140,  -130,  -210,  -250,   130,   -50,  -140,  -130],
    /* [UG] */ [  -60,  -100,  -140,  -150,   -50,    30,   -60,  -100],
    /* [AP] */ [ -280,  -162,  -277,  -329,  -140,   -60,  -280,  -162],
    /* [PA] */ [ -274,  -210,  -220,  -249,  -130,  -100,  -274,  -210],
];

pub static STACKPARAMS_ENTH: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA */
    /* [AU] */ [ -940,  -680, -1050, -1140,  -880,  -320, -2208, -2694],
    /* [UA] */ [ -680,  -770, -1040, -1240, -1280,  -700, -2081, -1247],
    /* [CG] */ [-1050, -1040, -1060, -1340, -1210,  -560, -1623, -1119],
    /* [GC] */ [-1140, -1240, -1340, -1490, -1260,  -830, -2407, -1729],
    /* [GU] */ [ -880, -1280, -1210, -1260, -1460, -1350,  -880, -1280],
    /* [UG] */ [ -700,  -320,  -830,  -560, -1350,  -930,  -320,  -700],
    /* [AP] */ [-2208, -2081, -1623, -2407,  -880,  -320, -2208, -2081],
    /* [PA] */ [-2694, -1247, -1119, -1729, -1280,  -700, -2694, -1247],
];

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK_EN37: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < E {
        let mut j = 0;
        while j < E {
            full[i][j] = Some(STACKPARAMS_EN37[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK_ENTH: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < E {
        let mut j = 0;
        while j < E {
            full[i][j] = Some(STACKPARAMS_ENTH[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};
