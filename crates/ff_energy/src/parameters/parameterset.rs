
pub const MAX_LOOP: usize = 30;
pub const B: usize = 4;  // A,C,G,U
pub const P: usize = 6;  // AU,UA,CG,GC,GU,UG
                         
pub type StackParams = [[i32; P]; P];
pub type LoopParams = [i32; MAX_LOOP + 1];
pub type MismatchParams = [[[i32; B]; B]; P];
pub type DangleParams = [[i32; B]; P];

pub type Int11Params = [[[[i32; B]; B]; P]; P];
pub type Int21Params = [[[[[i32; B]; B]; B]; P]; P];
pub type Int22Params = [[[[[[i32; B]; B]; B]; B]; P]; P];


pub struct ThermoParams {
    pub stack_en37: &'static StackParams,
    pub stack_enth: &'static StackParams,

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
    pub terminal_ru_en37: i32,
    pub terminal_ru_enth: i32,
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

    pub triloops: &'static [LoopEntry],
    pub tetraloops: &'static [LoopEntry],
    pub hexaloops: &'static [LoopEntry],
}

#[derive(Clone, Debug)]
pub struct LoopEntry {
    pub seq: &'static str,
    pub g37: i32,
    pub h: i32,
}

impl LoopEntry {
    #[inline]
    pub fn rescaled(&self, scale: f64) -> Self {
        let g37 = self.g37 as f64;
        let h = self.h as f64;
        Self {
            seq: self.seq,
            g37: (h - (h - g37) * scale).round() as i32,
            h: self.h,
        }
    }
}

