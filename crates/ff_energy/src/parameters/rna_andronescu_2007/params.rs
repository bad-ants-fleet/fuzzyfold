use crate::parameters::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_stacks::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_mimas::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_dangles::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int11::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int21::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_int22::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_loops::*;
use crate::parameters::rna_andronescu_2007::andronescu2007_hairpins::*;

pub struct AndronescuParams {
    pub stack: &'static ExtendedStackParams,
    pub mismatch_hairpin: &'static MismatchParams,
    pub mismatch_interior: &'static MismatchParams,
    pub mismatch_interior_1n: &'static MismatchParams,
    pub mismatch_interior_23: &'static MismatchParams,
    pub mismatch_multi: &'static MismatchParams,
    pub mismatch_exterior: &'static MismatchParams,
    pub dangle5: &'static DangleParams,
    pub dangle3: &'static DangleParams,
    pub int11: &'static Int11Params,
    pub int21: &'static Int21Params,
    pub int22: &'static Int22Params,
    pub hairpin: &'static LoopParams,
    pub bulge: &'static LoopParams,
    pub interior: &'static LoopParams,
    // Misc parameters
    pub duplex_init: i32,
    pub terminal_ru: i32,
    pub lxc: f64,
    // NINIO parameters
    pub ninio: i32,
    pub ninio_max: i32,
    // Multi-loop parameters
    pub ml_base: i32,
    pub ml_closing: i32,
    pub ml_intern: i32,
    // Special haipin parameters
    pub triloops: &'static [LoopEntry],
    pub tetraloops: &'static [LoopEntry],
    pub hexaloops: &'static [LoopEntry],
}

pub static RNA_ANDRONESCU_2007: AndronescuParams = AndronescuParams {
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


