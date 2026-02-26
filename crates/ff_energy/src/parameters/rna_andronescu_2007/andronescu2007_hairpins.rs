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

pub static TRILOOPS: &[LoopEntry] = &[ ];

pub static TETRALOOPS: &[LoopEntry] = &[
	loop_entry!(b"GGGGAC", 151), 
	loop_entry!(b"GGUGAC",  50),
	loop_entry!(b"CGAAAG", 147),
	loop_entry!(b"GGAGAC", 151),
	loop_entry!(b"CGCAAG", 147),
	loop_entry!(b"GGAAAC", 151),
	loop_entry!(b"CGGAAG", 147),
	loop_entry!(b"CUUCGG",  52),
	loop_entry!(b"CGUGAG", 147),
	loop_entry!(b"CGAAGG", 165),
	loop_entry!(b"CUACGG", 102),
	loop_entry!(b"GGCAAC",  60),
	loop_entry!(b"CGCGAG", 170),
	loop_entry!(b"UGAGAG",  75),
	loop_entry!(b"CGAGAG", 247),
	loop_entry!(b"AGAAAU", 202),
	loop_entry!(b"CGUAAG", 211),
	loop_entry!(b"CUAACG", 227),
	loop_entry!(b"UGAAAG", 278),
	loop_entry!(b"GGAAGC", 228),
	loop_entry!(b"GGGAAC", 301),
	loop_entry!(b"UGAAAA", 320),
	loop_entry!(b"AGCAAU", 272),
	loop_entry!(b"AGUAAU", 201),
	loop_entry!(b"CGGGAG", 261),
	loop_entry!(b"AGUGAU",  72),
	loop_entry!(b"GGCGAC", 101),
	loop_entry!(b"GGGAGC", 294),
	loop_entry!(b"GUGAAC", 358),
	loop_entry!(b"UGGAAA", 328),
];

pub static HEXALOOPS: &[LoopEntry] = &[ ];

