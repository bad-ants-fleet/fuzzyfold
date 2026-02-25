use crate::parameters::parameterset::LoopEntry;

pub static TRILOOPS_TURNER_2004: &[LoopEntry] = &[
    LoopEntry {seq: "CAACG", g37: 680, h: 2370},
    LoopEntry {seq: "CUUAC", g37: 690, h: 1080},
];

pub static TETRALOOPS_TURNER_2004: &[LoopEntry] = &[
    LoopEntry {seq: "CAACGG", g37: 550, h:   690}, 
    LoopEntry {seq: "CCAAGG", g37: 330, h: -1030},
    LoopEntry {seq: "CCACGG", g37: 370, h:  -330},
    LoopEntry {seq: "CCCAGG", g37: 340, h:  -890},
    LoopEntry {seq: "CCGAGG", g37: 350, h:  -660},
    LoopEntry {seq: "CCGCGG", g37: 360, h:  -750},
    LoopEntry {seq: "CCUAGG", g37: 370, h:  -350},
    LoopEntry {seq: "CCUCGG", g37: 250, h: -1390},
    LoopEntry {seq: "CUAAGG", g37: 360, h:  -760},
    LoopEntry {seq: "CUACGG", g37: 280, h: -1070},
    LoopEntry {seq: "CUCAGG", g37: 370, h:  -660},
    LoopEntry {seq: "CUCCGG", g37: 270, h: -1290},
    LoopEntry {seq: "CUGCGG", g37: 280, h: -1070},
    LoopEntry {seq: "CUUAGG", g37: 350, h:  -620},
    LoopEntry {seq: "CUUCGG", g37: 370, h: -1530},
    LoopEntry {seq: "CUUUGG", g37: 370, h:  -680},
];

pub static HEXALOOPS_TURNER_2004: &[LoopEntry] = &[
    LoopEntry { seq: "ACAGUACU", g37: 280, h: -1680},  
    LoopEntry { seq: "ACAGUGAU", g37: 360, h: -1140}, 
    LoopEntry { seq: "ACAGUGCU", g37: 290, h: -1280}, 
    LoopEntry { seq: "ACAGUGUU", g37: 180, h: -1540}, 
];

