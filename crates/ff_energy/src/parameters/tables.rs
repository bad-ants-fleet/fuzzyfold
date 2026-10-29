use crate::Base;

/// The temperature of the parameter sets, in degrees Celsius.
pub(crate) const T_REF_CELSIUS: f64 = 37.0;
/// 0 degrees Celsius in Kelvin.
pub(crate) const ZERO_CELSIUS_KELVIN: f64 = 273.15;
/// The longest tabulated loop; longer loops are extrapolated.
pub const MAX_LOOP: usize = 30;
/// Bases of the tables: A, C, G, U.
pub const B: usize = 4;
/// Canonical pairs of the tables: AU, UA, CG, GC, GU, UG.
pub const P: usize = 6;

/// Stacking energies, by the outer and the inner pair.
pub type StackParams = [[i32; P]; P];
/// Initiation energies by loop size (hairpin, bulge, interior).
pub type LoopParams = [i32; MAX_LOOP + 1];
/// Terminal mismatches, by pair and the two neighboring bases.
pub type MismatchParams = [[[i32; B]; B]; P];
/// Dangling ends, by pair and the neighboring base.
pub type DangleParams = [[i32; B]; P];

/// Interior loops of 1x1, 2x1 and 2x2 unpaired bases, by both pairs and
/// the unpaired bases.
pub type Int11Params = [[[[i32; B]; B]; P]; P];
pub type Int21Params = [[[[[i32; B]; B]; B]; P]; P];
pub type Int22Params = [[[[[[i32; B]; B]; B]; B]; P]; P];

/// Pairs of the extended tables: the canonical ones, and AP, PA, IU, UI,
/// CI, IC with pseudouridine (P) and inosine (I).
pub const E: usize = 12;
/// Stacking energies of the extended pairs (None: not parameterized).
pub type ExtendedStackParams = [[Option<i32>; E]; E];
/// Terminal penalties of the extended pairs.
pub type ClosingParams = [i32; E];

/// A special hairpin loop (its bases, closing pair included) and its
/// energy.
#[derive(Copy, Clone, Debug)]
pub struct LoopEntry {
    pub seq: &'static [Base],
    pub val: i32,
}

/// A thermodynamic parameter set: free energies at 37 C (`_en37`) and
/// enthalpies (`_enth`) for rescaling to other temperatures.
pub struct RNAThermoParams {
    pub stack_en37: &'static ExtendedStackParams,
    pub stack_enth: &'static ExtendedStackParams,
    pub mismatch_hairpin_en37: &'static MismatchParams,
    pub mismatch_hairpin_enth: &'static MismatchParams,
    pub mismatch_interior_en37: &'static MismatchParams,
    pub mismatch_interior_enth: &'static MismatchParams,
    pub mismatch_interior_1n_en37: &'static MismatchParams,
    pub mismatch_interior_1n_enth: &'static MismatchParams,
    pub mismatch_interior_23_en37: &'static MismatchParams,
    pub mismatch_interior_23_enth: &'static MismatchParams,
    pub mismatch_multi_en37: &'static MismatchParams,
    pub mismatch_multi_enth: &'static MismatchParams,
    pub mismatch_exterior_en37: &'static MismatchParams,
    pub mismatch_exterior_enth: &'static MismatchParams,
    pub dangle5_en37: &'static DangleParams,
    pub dangle5_enth: &'static DangleParams,
    pub dangle3_en37: &'static DangleParams,
    pub dangle3_enth: &'static DangleParams,
    pub int11_en37: &'static Int11Params,
    pub int11_enth: &'static Int11Params,
    pub int21_en37: &'static Int21Params,
    pub int21_enth: &'static Int21Params,
    pub int22_en37: &'static Int22Params,
    pub int22_enth: &'static Int22Params,
    pub hairpin_en37: &'static LoopParams,
    pub hairpin_enth: &'static LoopParams,
    pub bulge_en37: &'static LoopParams,
    pub bulge_enth: &'static LoopParams,
    pub interior_en37: &'static LoopParams,
    pub interior_enth: &'static LoopParams,
    // Misc parameters
    pub duplex_init_en37: i32,
    pub duplex_init_enth: i32,
    pub closing_en37: &'static ClosingParams,
    pub closing_enth: &'static ClosingParams,
    pub lxc: f64,
    // NINIO parameters
    pub ninio_en37: i32,
    pub ninio_enth: i32,
    pub ninio_max: i32,
    // Multi-loop parameters
    pub ml_base_en37: i32,
    pub ml_base_enth: i32,
    pub ml_closing_en37: i32,
    pub ml_closing_enth: i32,
    pub ml_intern_en37: i32,
    pub ml_intern_enth: i32,
    // Special haipin parameters
    pub triloops_en37: &'static [LoopEntry],
    pub triloops_enth: &'static [LoopEntry],
    pub tetraloops_en37: &'static [LoopEntry],
    pub tetraloops_enth: &'static [LoopEntry],
    pub hexaloops_en37: &'static [LoopEntry],
    pub hexaloops_enth: &'static [LoopEntry],
}



