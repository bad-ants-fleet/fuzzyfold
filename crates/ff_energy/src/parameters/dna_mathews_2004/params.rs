use crate::parameters::ThermoParams;
use crate::parameters::dna_mathews_2004::mathews2004_stacks::*;
use crate::parameters::dna_mathews_2004::mathews2004_mimas::*;
use crate::parameters::dna_mathews_2004::mathews2004_int11::*;
use crate::parameters::dna_mathews_2004::mathews2004_int21::*;
use crate::parameters::dna_mathews_2004::mathews2004_int22::*;
use crate::parameters::dna_mathews_2004::mathews2004_loops::*;
use crate::parameters::dna_mathews_2004::mathews2004_dangles::*;
use crate::parameters::parameterset::LoopEntry;

pub static LOOPS: &[LoopEntry] = &[ ];

pub static DNA_MATHEWS_2004: ThermoParams = ThermoParams {
    stack_en37: &STACK_EN37,
    stack_enth: &STACK_ENTH,

    mismatch_hairpin_en37: &MISMATCH_HAIRPIN_EN37,
    mismatch_hairpin_enth: &MISMATCH_HAIRPIN_ENTH,
    mismatch_interior_en37: &MISMATCH_INTERIOR_EN37,
    mismatch_interior_enth: &MISMATCH_INTERIOR_ENTH,
    mismatch_interior_1n_en37: &MISMATCH_INTERIOR_1N_EN37,
    mismatch_interior_1n_enth: &MISMATCH_INTERIOR_1N_ENTH,
    mismatch_interior_23_en37: &MISMATCH_INTERIOR_23_EN37,
    mismatch_interior_23_enth: &MISMATCH_INTERIOR_23_ENTH,
    mismatch_multi_en37: &MISMATCH_MULTI_EN37,
    mismatch_multi_enth: &MISMATCH_MULTI_ENTH,
    mismatch_exterior_en37: &MISMATCH_EXTERIOR_EN37,
    mismatch_exterior_enth: &MISMATCH_EXTERIOR_ENTH,
                                                          
    dangle5_en37: &DANGLE5_EN37,
    dangle5_enth: &DANGLE5_ENTH,
    dangle3_en37: &DANGLE3_EN37,
    dangle3_enth: &DANGLE3_ENTH,
                                                                              
    int11_en37: &INT11_EN37,
    int11_enth: &INT11_ENTH,
    int21_en37: &INT21_EN37,
    int21_enth: &INT21_ENTH,
    int22_en37: &INT22_EN37,
    int22_enth: &INT22_ENTH,
             
    hairpin_en37: &HAIRPIN_EN37,
    hairpin_enth: &HAIRPIN_ENTH,
    bulge_en37: &BULGE_EN37,
    bulge_enth: &BULGE_ENTH,
    interior_en37: &INTERIOR_EN37,
    interior_enth: &INTERIOR_ENTH,

    duplex_init_en37:  100,
    duplex_init_enth: -720,
    terminal_ru_en37:    0,
    terminal_ru_enth:  320,
    lxc: 107.856,

    ninio_en37: 40,
    ninio_enth:  0,
    ninio_max: 300,

    ml_base_en37: 20,
    ml_base_enth:  0,
    ml_closing_en37: 300,
    ml_closing_enth: 900,
    ml_intern_en37:   20,
    ml_intern_enth:    0,

    triloops_en37: LOOPS,
    triloops_enth: LOOPS,
    tetraloops_en37: LOOPS,
    tetraloops_enth: LOOPS,
    hexaloops_en37: LOOPS,
    hexaloops_enth: LOOPS,
};


