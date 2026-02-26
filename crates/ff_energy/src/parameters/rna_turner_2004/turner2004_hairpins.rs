use crate::Base;
use crate::parameters::parameterset::LoopEntry;

const fn convert<const N: usize>(bytes: &[u8; N]) -> [Base; N] {
    let mut out = [Base::A; N];
    let mut i = 0;
    while i < N {
        out[i] = 
            match bytes[i] {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'U' => Base::U,
                _ => panic!("Invalid base"),
            };
        i += 1;
    }
    out
}

macro_rules! loop_entry {
    ($seq:literal, $val:expr) => {
        LoopEntry {
            seq: &convert($seq),
            val: $val,
        }
    };
}

pub static TRILOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"CAACG", 680),
    loop_entry!(b"CUUAC", 690),
];

pub static TRILOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"CAACG", 2370),
    loop_entry!(b"CUUAC", 1080),
];

pub static TETRALOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"CAACGG", 550), 
    loop_entry!(b"CCAAGG", 330),
    loop_entry!(b"CCACGG", 370),
    loop_entry!(b"CCCAGG", 340),
    loop_entry!(b"CCGAGG", 350),
    loop_entry!(b"CCGCGG", 360),
    loop_entry!(b"CCUAGG", 370),
    loop_entry!(b"CCUCGG", 250),
    loop_entry!(b"CUAAGG", 360),
    loop_entry!(b"CUACGG", 280),
    loop_entry!(b"CUCAGG", 370),
    loop_entry!(b"CUCCGG", 270),
    loop_entry!(b"CUGCGG", 280),
    loop_entry!(b"CUUAGG", 350),
    loop_entry!(b"CUUCGG", 370),
    loop_entry!(b"CUUUGG", 370),
];

pub static TETRALOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"CAACGG",   690), 
    loop_entry!(b"CCAAGG", -1030),
    loop_entry!(b"CCACGG",  -330),
    loop_entry!(b"CCCAGG",  -890),
    loop_entry!(b"CCGAGG",  -660),
    loop_entry!(b"CCGCGG",  -750),
    loop_entry!(b"CCUAGG",  -350),
    loop_entry!(b"CCUCGG", -1390),
    loop_entry!(b"CUAAGG",  -760),
    loop_entry!(b"CUACGG", -1070),
    loop_entry!(b"CUCAGG",  -660),
    loop_entry!(b"CUCCGG", -1290),
    loop_entry!(b"CUGCGG", -1070),
    loop_entry!(b"CUUAGG",  -620),
    loop_entry!(b"CUUCGG", -1530),
    loop_entry!(b"CUUUGG",  -680),
];

pub static HEXALOOPS_EN37: &[LoopEntry] = &[
    loop_entry!(b"ACAGUACU", 280),  
    loop_entry!(b"ACAGUGAU", 360), 
    loop_entry!(b"ACAGUGCU", 290), 
    loop_entry!(b"ACAGUGUU", 180), 
];

pub static HEXALOOPS_ENTH: &[LoopEntry] = &[
    loop_entry!(b"ACAGUACU", -1680),  
    loop_entry!(b"ACAGUGAU", -1140), 
    loop_entry!(b"ACAGUGCU", -1280), 
    loop_entry!(b"ACAGUGUU", -1540), 
];

