use crate::parameters::RNAThermoParams;
use crate::parameters::rna_turner_2004::turner2004_mimas::*;
use crate::parameters::rna_turner_2004::turner2004_int11::*;
use crate::parameters::rna_turner_2004::turner2004_int21::*;
use crate::parameters::rna_turner_2004::turner2004_int22::*;
use crate::parameters::rna_turner_2004::turner2004_loops::*;
use crate::parameters::rna_turner_2004::turner2004_dangles::*;
use crate::parameters::rna::extended_stacks::*;
use crate::parameters::rna::extended_hairpins::*;

pub static RNA_EXTENDED: RNAThermoParams = RNAThermoParams {
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

    duplex_init_en37:  410,
    duplex_init_enth:  360,
    terminal_ru_en37:   50,
    terminal_ru_enth:  370,
    terminal_ap_en37:   31,
    terminal_ap_enth: -204,
    lxc: 107.856,

    ninio_en37:  60,
    ninio_enth: 320,
    ninio_max:  300,

    ml_base_en37: 0,
    ml_base_enth: 0,
    ml_closing_en37:  930,
    ml_closing_enth: 3000,
    ml_intern_en37:  -90,
    ml_intern_enth: -220,

    triloops_en37: COLLECTED_TRILOOPS_EN37,
    triloops_enth: COLLECTED_TRILOOPS_ENTH,
    tetraloops_en37: COLLECTED_TETRALOOPS_EN37,
    tetraloops_enth: COLLECTED_TETRALOOPS_ENTH,
    hexaloops_en37: COLLECTED_HEXALOOPS_EN37,
    hexaloops_enth: COLLECTED_HEXALOOPS_ENTH,
};


