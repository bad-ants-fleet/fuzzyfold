use crate::parameters::parameterset::StackParams;
use crate::parameters::parameterset::ExtendedStackParams;
use crate::parameters::parameterset::{P, E};

pub static STACKPARAMS: StackParams = [
    /* [cl] [i]:      AU     UA     CG     GC     GU     UG */
    /* [AU] */ [  -74,   -71,  -154,  -132,   -91,   -12],
    /* [UA] */ [  -71,   -77,  -149,  -139,   -69,   -12],
    /* [CG] */ [ -154,  -149,  -152,  -216,  -168,   -58],
    /* [GC] */ [ -132,  -139,  -216,  -216,  -157,   -97],
    /* [GU] */ [  -91,   -69,  -168,  -157,   -29,   -73],
    /* [UG] */ [  -12,   -12,   -58,   -97,   -73,   -58],
];

/// Full ExtendedStackParams built from STACKPARAMS.
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
