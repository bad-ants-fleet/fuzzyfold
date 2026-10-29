//! # Pseudouridine Parameters
//!
//! **Authors:** Graham A. Hudson, Richard J. Bloomingdale, and Brent M. Znosko  
//! **Title:** *Thermodynamic contribution and nearest-neighbor parameters of
//! pseudouridine-adenosine base pairs in oligoribonucleotides*  
//! **Journal:** RNA 19:1474-1482  
//! **Year:** 2013  
//! **DOI:** 10.1261/rna.039610.113
//! 
//! # Inosine Parameters 
//!   
//! For `[IU][IU]`; `[IU][UI]`; `[UI][IU]`; `[UI][UI]`:
//! 
//! **Authors:** Daniel J. Wright, Jamie L. Rice, Dawn M. Yanker, and Brent M. Znosko
//! **Title:** Nearest Neighbor Parameters for Inosine-Uridine Pairs in RNA Duplexes
//! **Journal:** Biochemistry, 46, 4625-4634
//! **Year:** 2007 
//! **DOI:** 10.1021/bi0616910
//!  
//! **Authors:** Daniel J. Wright, Christopher R. Force and Brent M. Znosko
//! **Title:** Stability of RNA duplexes containing inosine-cytosine pairs
//! **Journal:** Nucleic Acids Research, Vol. 46, No. 22 12099-12108
//! **Year:** 2018
//! **DOI:** 10.1093/nar/gky907
//! 
//! For the rest: 
//! 
//! **Authors:** Hamid Reza Mansouri Khosravi, Thomas Spicher, Cornelia Vesely, Ronny Lorenz, Karolina Bartosik, Julia Thaler, Victorio Jauregui-Matos, Ronald Micura, Ivo L. Hofacker, Michael F. Jantsch
//! **Title:** Systematic measurement of thermodynamic nearest-neighbor parameters for Inosine-containing double-stranded RNAs using a fluorophore-quencher-based approach
//! **Journal:** NAR, submitted
//! **Year:** 2026
//! **DOI:** 
//!
//! ## Description
//!
//! This module provides stacking parameters including **P** for
//! pseudouridine and **I** for inosine.
//!
//! ## Implementation Notes
//!
//! - All unspecified pseudouridine interactions are treated as **U**.
//! - `[AP][AP]` = `[AP][AU]`
//! - `[AP][GU]` = `[AU][GU]`
//! - `[GP]` = `[GU]` (not even shown in table.)
//! - All unspecified inosine interactions are treated as **G**.
//! - `[IU][IU]` = `[IU][GU]`
//! - `[AP][IU]` = `[AP][GU]`
//! - `[IU][AP]` = `[IU][AU]`
//! 
//! 
use crate::parameters::tables::ExtendedStackParams;
use crate::parameters::tables::E;

pub static STACKPARAMS_EN37: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
    /* [AU] */ [ -110,   -90,  -210,  -220,  -140,   -60,  -280,  -274,   -65,    36,  -147,  -151],                  
    /* [UA] */ [  -90,  -130,  -210,  -240,  -130,  -100,  -162,  -210,   -77,     8,  -125,  -129],
    /* [CG] */ [ -210,  -210,  -240,  -330,  -210,  -140,  -277,  -220,  -146,   -90,  -135,  -236],
    /* [GC] */ [ -220,  -240,  -330,  -340,  -250,  -150,  -329,  -249,  -242,   -62,  -245,  -246],
    /* [GU] */ [ -140,  -130,  -210,  -250,   130,   -50,  -140,  -130,   163,    54,   -96,   -99],
    /* [UG] */ [  -60,  -100,  -140,  -150,   -50,    30,   -60,  -100,   -20,   146,   -15,   -26],
    /* [AP] */ [ -280,  -162,  -277,  -329,  -140,   -60,  -280,  -162,  -140,   -60,  -277,  -329],
    /* [PA] */ [ -274,  -210,  -220,  -249,  -130,  -100,  -274,  -210,  -130,  -100,  -220,  -249],
    /* [IU] */ [  -65,   -77,  -146,  -242,   163,   -20,   -65,   -77,   358,   266,  -146,  -242],           
    /* [UI] */ [   36,     8,   -90,   -62,    54,   146,    36,     8,   266,   223,   -90,   -62],
    /* [CI] */ [ -147,  -125,  -135,  -245,   -96,   -15,  -147,  -125,   -96,   -15,  -135,  -245],           
    /* [IC] */ [ -151,  -129,  -236,  -246,   -99,   -26,  -151,  -129,   -99,   -26,  -236,  -246],
];

pub static STACKPARAMS_ENTH: [[i32; E]; E] = [
    /* [cl] [ri]:  AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC */
    /* [AU] */ [ -940,  -680, -1050, -1140,  -880,  -320, -2208, -2694, -1058,  -115, -1098, -1195],
    /* [UA] */ [ -680,  -770, -1040, -1240, -1280,  -700, -2081, -1247, -1274,  -572, -1148, -1096],
    /* [CG] */ [-1050, -1040, -1060, -1340, -1210,  -560, -1623, -1119,  -992,  -430,  -951, -1330],
    /* [GC] */ [-1140, -1240, -1340, -1490, -1260,  -830, -2407, -1729, -1620,  -431, -1345, -1333],
    /* [GU] */ [ -880, -1280, -1210, -1260, -1460, -1350,  -880, -1280, -1355,  -945, -1156, -1134],
    /* [UG] */ [ -320,  -700,  -560,  -830, -1350,  -930,  -320,  -700, -1265,  -619,  -391,  -569],
    /* [AP] */ [-2208, -2081, -1623, -2407,  -880,  -320, -2208, -2081,  -880,  -320, -1623, -2470],
    /* [PA] */ [-2694, -1247, -1119, -1729, -1280,  -700, -2694, -1247, -1280,  -700, -1119, -1729],
    /* [IU] */ [-1058, -1274,  -992, -1620, -1355, -1265, -1058, -1274,   170,   953,  -992, -1620],           
    /* [UI] */ [ -115,  -572,  -430,  -431,  -945,  -619,  -115,  -572,   953,   841,  -430,  -431],
    /* [CI] */ [-1098, -1148,  -951, -1345, -1156,  -391, -1098, -1148, -1156,  -391,  -951, -1345],           
    /* [IC] */ [-1195, -1096, -1330, -1333, -1134,  -569, -1195, -1096, -1134,  -569, -1330, -1333],
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