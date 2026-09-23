use crate::parameters::rna_mt09::dangles::*;
use crate::parameters::rna_mt09::hairpins::*;
use crate::parameters::rna_mt09::int11::*;
use crate::parameters::rna_mt09::int21::*;
use crate::parameters::rna_mt09::int22::*;
use crate::parameters::rna_mt09::loops::*;
use crate::parameters::rna_mt09::mimas::*;
use crate::parameters::rna_mt09::stacks::*;
use crate::parameters::*;

pub static RNA_MT09: AndronescuParams = AndronescuParams {
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

    duplex_init: -46,
    terminal_ru: 50,
    lxc: 107.9,

    ninio: 50,
    ninio_max: 300,

    ml_base: 2,
    ml_closing: 339,
    ml_intern: 3,

    triloops: TETRALOOPS,
    tetraloops: TETRALOOPS,
    hexaloops: TETRALOOPS,
};
