use crate::parameters::tables::ClosingParams;

// AP, PA, IU, UI, CI, IC parameters are dummy parameters

pub static CLOSING_PEN37: ClosingParams = 
    /* [cl]:      AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [   0,      0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0];


pub static CLOSING_ENTH: ClosingParams =
    /* [cl]:       AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [  320,   320,     0,     0,   320,   320,   320,   320,   320,   320,     0,     0];