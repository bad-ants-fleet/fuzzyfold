use crate::parameters::parameterset::StackParams;
use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::{P, E};

pub static STACKPARAMS_EN37: StackParams = [
    /* [cl] [i]:   AU     UA     CG     GC     GU     UG */
    /* [AU] */ [ -110,   -90,  -210,  -220,  -140,   -60],
    /* [UA] */ [  -90,  -130,  -210,  -240,  -130,  -100],
    /* [CG] */ [ -210,  -210,  -240,  -330,  -210,  -140],
    /* [GC] */ [ -220,  -240,  -330,  -340,  -250,  -150],
    /* [GU] */ [ -140,  -130,  -210,  -250,   130,   -50],
    /* [UG] */ [  -60,  -100,  -140,  -150,   -50,    30],
];

pub static STACKPARAMS_ENTH: StackParams = [
    /* [cl] [i]:   AU     UA     CG     GC     GU     UG */
    /* [AU] */ [ -940,  -680, -1050, -1140,  -880,  -320],
    /* [UA] */ [ -680,  -770, -1040, -1240, -1280,  -700],
    /* [CG] */ [-1050, -1040, -1060, -1340, -1210,  -560],
    /* [GC] */ [-1140, -1240, -1340, -1490, -1260,  -830],
    /* [GU] */ [ -880, -1280, -1210, -1260, -1460, -1350],
    /* [UG] */ [ -320,  -700,  -560,  -830, -1350,  -930],
];

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK_EN37: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < P {
        let mut j = 0;
        while j < P {
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
    while i < P {
        let mut j = 0;
        while j < P {
            full[i][j] = Some(STACKPARAMS_ENTH[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};
