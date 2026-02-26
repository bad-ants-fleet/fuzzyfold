use crate::parameters::parameterset::MismatchParams;

pub static MISMATCH_HAIRPIN_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_HAIRPIN_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_1N_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_1N_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_23_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_INTERIOR_23_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_MULTI_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_MULTI_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_EXTERIOR_EN37: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

pub static MISMATCH_EXTERIOR_ENTH: MismatchParams = [
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [AU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [AU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UA][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UA][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [CG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [CG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GC][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GC][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [GU][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [GU][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
    [ /* [cl][5] [3]:    A      C      G      U */
      /* [UG][A] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][C] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][G] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
      /* [UG][U] */ [i32::MAX, i32::MAX, i32::MAX, i32::MAX],
    ],
];

