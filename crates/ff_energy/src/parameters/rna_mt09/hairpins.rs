use crate::Base;
use crate::parameters::parameterset::LoopEntry;

const fn convert<const N: usize>(bytes: &[u8; N]) -> [Base; N] {
    let mut out = [Base::A; N];
    let mut i = 0;
    while i < N {
        out[i] = match bytes[i] {
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
        LoopEntry { seq: &convert($seq), val: $val }
    };
}

pub static TRILOOPS: &[LoopEntry] = &[];

pub static TETRALOOPS: &[LoopEntry] = &[
	loop_entry!(b"GGGGAC", 150),
	loop_entry!(b"GGUGAC", 12),
	loop_entry!(b"CGAAAG", 85),
	loop_entry!(b"GGAGAC", 153),
	loop_entry!(b"CGCAAG", 49),
	loop_entry!(b"GGAAAC", 47),
	loop_entry!(b"CGGAAG", 57),
	loop_entry!(b"CUUCGG", 4),
	loop_entry!(b"CGUGAG", 19),
	loop_entry!(b"CGAAGG", 46),
	loop_entry!(b"CUACGG", 2),
	loop_entry!(b"GGCAAC", 76),
	loop_entry!(b"CGCGAG", 73),
	loop_entry!(b"UGAGAG", 40),
	loop_entry!(b"CGAGAG", 235),
	loop_entry!(b"AGAAAU", 118),
	loop_entry!(b"CGUAAG", 113),
	loop_entry!(b"CUAACG", 91),
	loop_entry!(b"UGAAAG", 159),
	loop_entry!(b"GGAAGC", 138),
	loop_entry!(b"GGGAAC", 135),
	loop_entry!(b"UGAAAA", 172),
	loop_entry!(b"AGCAAU", 169),
	loop_entry!(b"AGUAAU", 201),
	loop_entry!(b"CGGGAG", 146),
	loop_entry!(b"AGUGAU", 185),
	loop_entry!(b"GGCGAC", 197),
	loop_entry!(b"GGGAGC", 192),
	loop_entry!(b"GUGAAC", 257),
	loop_entry!(b"UGGAAA", 245),
];

pub static HEXALOOPS: &[LoopEntry] = &[];

