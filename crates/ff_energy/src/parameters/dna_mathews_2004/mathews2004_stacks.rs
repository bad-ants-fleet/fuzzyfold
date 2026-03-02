use crate::parameters::parameterset::StackParams;
use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::{P, E};

pub static STACKPARAMS_EN37: StackParams = [
    /* [cl] [i]:   AT     TA     CG     GC     GT     TG */
    /* [AT] */ [  -90,  -100,  -130,  -140,    10,    70],
    /* [TA] */ [ -100,   -60,  -150,  -130,    30,    40],
    /* [CG] */ [ -130,  -150,  -220,  -180,   -30,   -50],
    /* [GC] */ [ -140,  -130,  -180,  -220,   -60,    10],
    /* [GT] */ [   10,    30,   -30,   -60,   120,    70],
    /* [TG] */ [   70,    40,   -50,    10,    70,    50],
];

pub static STACKPARAMS_ENTH: StackParams = [
    /* [cl] [i]:   AT     TA     CG     GC     GT     TG */
    /* [AT] */ [ -530,  -720,  -580,  -780,  -170,   370],
    /* [TA] */ [ -720,  -670,  -990,  -850,   -90,   -90],
    /* [CG] */ [ -580,  -990,  -980,  -750,  -240,  -430],
    /* [GC] */ [ -780,  -850,  -750,  -790,  -240,   430],
    /* [GT] */ [ -170,   -90,  -240,  -240,   500,   800],
    /* [TG] */ [  340,   -90,  -430,   430,   800,  -670],
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
