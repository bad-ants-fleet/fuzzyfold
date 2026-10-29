use crate::parameters::tables::ClosingParams;

pub static CLOSING_PEN37: ClosingParams = 
    /* [cl]:      AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [  50,    50,     0,     0,    50,    50,    31,    31,   -21,   -21,   -56,   -56];


pub static CLOSING_ENTH: ClosingParams =
    /* [cl]:       AU     UA     CG     GC     GU     UG     AP     PA     IU     UI     CI     IC*/
               [  370,   370,     0,     0,   370,   370,  -204,  -204,   -34,   -34,  -255,  -255];
