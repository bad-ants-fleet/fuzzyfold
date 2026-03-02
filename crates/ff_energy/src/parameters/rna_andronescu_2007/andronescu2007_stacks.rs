use crate::parameters::parameterset::StackParams;
use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::{P, E};

pub static STACKPARAMS: StackParams = [
    /* [cl] [i]:   AU     UA     CG     GC     GU     UG */
    /* [AU] */ [  -99,   -69,  -199,  -189,   -88,    -1],
    /* [UA] */ [  -69,   -91,  -178,  -211,   -47,     0],
    /* [CG] */ [ -199,  -178,  -203,  -271,  -178,   -85],
    /* [GC] */ [ -189,  -211,  -271,  -300,  -193,  -127],
    /* [GU] */ [  -88,   -47,  -178,  -193,    30,   -71],
    /* [UG] */ [   -1,     0,   -85,  -127,   -71,   -70],
];

/// The parameters embedded into whatever dimension is currently 
/// used by the ViennaRNA-style stacking tables.
pub const STACK: ExtendedStackParams = {
    let mut full: ExtendedStackParams = [[None; E]; E];

    let mut i = 0;
    while i < P {
        let mut j = 0;
        while j < P {
            full[i][j] = Some(STACKPARAMS[i][j]);
            j += 1;
        }
        i += 1;
    }

    full
};
