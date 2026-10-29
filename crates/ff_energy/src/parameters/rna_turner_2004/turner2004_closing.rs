use crate::parameters::tables::ClosingParams;

// AP, PA, IU, UI, CI, IC parameters are dummy parameters


pub static CLOSING_PEN37: ClosingParams = 
    /* [cl]:      AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [  50,    50,     0,     0,    50,    50,    50,    50,    50,    50,     0,     0];


pub static CLOSING_ENTH: ClosingParams =
    /* [cl]:       AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [  370,   370,     0,     0,   370,   370,   370,   370,   370,   370,     0,     0];


   