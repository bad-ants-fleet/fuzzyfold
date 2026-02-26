use crate::parameters::FittedParams;
use crate::parameters::rna_andronescu_2007::andronescu2007_stacks::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_mimas::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_dangles::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int11::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int21::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int22::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_loops::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_hairpins::*;

pub static RNA_ANDRONESCU_2007: FittedParams = FittedParams {
    stack: &STACK,
    mismatch_hairpin: &MISMATCH_HAIRPIN,
    mismatch_interior: &MISMATCH_INTERIOR,
    mismatch_interior_1n: &MISMATCH_INTERIOR_1N,
    mismatch_interior_23: &MISMATCH_INTERIOR_23,
    mismatch_multi: &MISMATCH_MULTI,
    mismatch_exterior: &MISMATCH_EXTERIOR,
                                                          
    dangle5: &DANGLE5,
    dangle3: &DANGLE3,
                                                                              
    int11: &INT11,
    int21: &INT21,
    int22: &INT22,
             
    hairpin: &HAIRPIN,
    bulge: &BULGE,
    interior: &INTERIOR,

    duplex_init: 410,
    terminal_ru:  11,
    lxc: 107.856,

    ninio:  50,
    ninio_max:  300,

    ml_base: 4,
    ml_closing:  440,
    ml_intern:  3,

    triloops: TETRALOOPS,
    tetraloops: TETRALOOPS,
    hexaloops: TETRALOOPS,
};


