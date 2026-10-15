use crate::parameters::parameterset::DangleParams;

pub static DANGLE5: DangleParams = [
    /* [cl] [D]:    A      C      G      U */
    /* [AU] */ [  -32,   -70,   -90,   -35],
    /* [UA] */ [  -32,   -42,   -68,   -45],
    /* [CG] */ [  -87,   -24,   -98,   -86],
    /* [GC] */ [  -46,   -37,   -95,   -54],
    /* [GU] */ [  -32,   -87,  -108,   -11],
    /* [UG] */ [  -32,   -19,   -86,   -11],
];

pub static DANGLE3: DangleParams = [
    /* [cl] [D]:    A      C      G      U */
    /* [AU] */ [  -14,     0,     0,    -7],
    /* [UA] */ [  -32,    -7,    -6,    -6],
    /* [CG] */ [  -32,   -19,   -13,    -5],
    /* [GC] */ [   -9,   -17,     0,    -1],
    /* [GU] */ [    0,    -1,     0,   -11],
    /* [UG] */ [    0,     0,     0,     0],
];

