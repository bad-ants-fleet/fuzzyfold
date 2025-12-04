//! A parser for the parameter file format shipped with ViennaRNA.
//!
//! IMPORTANT: This module provides the hardcoded indices of parameter entries
//! as specified in the parsed file! Verify that orders are correct, before
//! using this parser!!
//!
//! No error handling implemented. If you got a malformed parameter file, 
//! it's going to panic through .expect()
//!

/// Importing only EnergyTables and it's indices.
use crate::{Base, NucleotideVec};
use crate::PairTypeRNA;
use crate::nn_parsing::EnergyTables;
use itertools::Itertools;

const PARAM_FILE_PAIR_ORDER: [PairTypeRNA; 7] = [
    PairTypeRNA::CG,
    PairTypeRNA::GC,
    PairTypeRNA::GU,
    PairTypeRNA::UG,
    PairTypeRNA::AU,
    PairTypeRNA::UA,
    PairTypeRNA::NN,
];

const PARAM_FILE_MM_ORDER: [Base; 5] = [
    Base::N,
    Base::A,
    Base::C,
    Base::G,
    Base::U,
];

pub trait SectionParser {
    fn parse_line(&mut self, line: &str, tables: &mut EnergyTables);
}

macro_rules! impl_stack_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            outer: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (inner, token) in line
                    .split_whitespace()
                    .take(PARAM_FILE_PAIR_ORDER.len())
                    .enumerate()
                {
                    let val = token.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            token
                    ));
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_PAIR_ORDER[inner] as usize;
                    tables.$field[i1][i2] = Some(val);
                }
                self.outer += 1;
            }
        }
    };
}

impl_stack_parser!(Stack, stack);
impl_stack_parser!(StackEnthalpies, stack_enthalpies);

macro_rules! impl_mismatch_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            outer: usize,
            m5: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (m3, token) in line.split_whitespace()
                    .take(PARAM_FILE_MM_ORDER.len()).enumerate() 
                {
                    let val = token.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            token
                    ));
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_MM_ORDER[self.m5] as usize;
                    let i3 = PARAM_FILE_MM_ORDER[m3] as usize;
                    tables.$field[i1][i2][i3] = Some(val);
                }
                self.m5 += 1;
                if self.m5 == PARAM_FILE_MM_ORDER.len() {
                    self.outer += 1;
                    self.m5 = 0;
                }
            }
        }
    };
}

impl_mismatch_parser!(MismatchHairpin, mismatch_hairpin);
impl_mismatch_parser!(MismatchHairpinEnthalpies, mismatch_hairpin_enthalpies);
impl_mismatch_parser!(MismatchInterior, mismatch_interior);
impl_mismatch_parser!(MismatchInteriorEnthalpies, mismatch_interior_enthalpies);
impl_mismatch_parser!(MismatchInterior1n, mismatch_interior_1n);
impl_mismatch_parser!(MismatchInterior1nEnthalpies, mismatch_interior_1n_enthalpies);
impl_mismatch_parser!(MismatchInterior23, mismatch_interior_23);
impl_mismatch_parser!(MismatchInterior23Enthalpies, mismatch_interior_23_enthalpies);
impl_mismatch_parser!(MismatchMulti, mismatch_multi);
impl_mismatch_parser!(MismatchMultiEnthalpies, mismatch_multi_enthalpies);
impl_mismatch_parser!(MismatchExterior, mismatch_exterior);
impl_mismatch_parser!(MismatchExteriorEnthalpies, mismatch_exterior_enthalpies);

macro_rules! impl_dangle_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            outer: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (m5, token) in line.split_whitespace()
                    .take(PARAM_FILE_MM_ORDER.len()).enumerate() 
                {
                    let val = if token == "INF" {
                        i32::MAX
                    } else {
                        token.parse::<i32>().unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                                stringify!($field),
                                line,
                                token
                            )
                        })
                    };
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_MM_ORDER[m5] as usize;
                    tables.$field[i1][i2] = Some(val);
                }
                self.outer += 1;
            }
        }
    };
}

impl_dangle_parser!(Dangle5, dangle5);
impl_dangle_parser!(Dangle5Enthalpies, dangle5_enthalpies);
impl_dangle_parser!(Dangle3, dangle3);
impl_dangle_parser!(Dangle3Enthalpies, dangle3_enthalpies);

macro_rules! impl_int11_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            outer: usize,
            inner: usize,
            mm5: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (mm3, token) in line.split_whitespace()
                    .take(PARAM_FILE_MM_ORDER.len()).enumerate() 
                {
                    let val = token.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            token
                    ));
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_PAIR_ORDER[self.inner] as usize;
                    let i3 = PARAM_FILE_MM_ORDER[self.mm5] as usize;
                    let i4 = PARAM_FILE_MM_ORDER[mm3] as usize;
                    tables.$field[i1][i2][i3][i4] = Some(val);
                }

                self.mm5 += 1;
                if self.mm5 == PARAM_FILE_MM_ORDER.len() {
                    self.mm5 = 0;
                    self.inner += 1;
                } 
                if self.inner == PARAM_FILE_PAIR_ORDER.len() {
                    self.outer += 1;
                    self.inner = 0;
                }
            }
        }
    };
}

impl_int11_parser!(Int11, int11);
impl_int11_parser!(Int11Enthalpies, int11_enthalpies);

macro_rules! impl_int21_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            outer: usize,
            inner: usize,
            mm55: usize,
            mm53: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (mm3, token) in line.split_whitespace()
                    .take(PARAM_FILE_MM_ORDER.len()).enumerate() 
                {
                    let val = token.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            token
                    ));
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_PAIR_ORDER[self.inner] as usize;
                    let i3 = PARAM_FILE_MM_ORDER[self.mm55] as usize;
                    let i4 = PARAM_FILE_MM_ORDER[self.mm53] as usize;
                    let i5 = PARAM_FILE_MM_ORDER[mm3] as usize;
                    tables.$field[i1][i2][i3][i4][i5] = Some(val);
                }
                self.mm53 += 1;
                if self.mm53 == PARAM_FILE_MM_ORDER.len() {
                    self.mm55 += 1;
                    self.mm53 = 0;
                }
                if self.mm55 == PARAM_FILE_MM_ORDER.len() {
                    self.mm55 = 0;
                    self.inner += 1;
                } 
                if self.inner == PARAM_FILE_PAIR_ORDER.len() {
                    self.outer += 1;
                    self.inner = 0;
                }
            }
        }
    };
}

impl_int21_parser!(Int21, int21);
impl_int21_parser!(Int21Enthalpies, int21_enthalpies);

macro_rules! impl_int22_parser {
    ($struct_name:ident, $field:ident) => {
        pub struct $struct_name {
            outer: usize,
            inner: usize,
            mm55: usize,
            mm53: usize,
            mm35: usize,
        }

        impl Default for $struct_name {
            fn default() -> Self {
                Self {
                    outer: 0,
                    inner: 0,
                    mm55: 1,
                    mm53: 1,
                    mm35: 1,
                }
            }
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                for (mm33, token) in line.split_whitespace()
                    .take(PARAM_FILE_MM_ORDER.len() - 1).enumerate() 
                {
                    let val = token.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            token
                    ));
                    let i1 = PARAM_FILE_PAIR_ORDER[self.outer] as usize;
                    let i2 = PARAM_FILE_PAIR_ORDER[self.inner] as usize;
                    let i3 = PARAM_FILE_MM_ORDER[self.mm55] as usize;
                    let i4 = PARAM_FILE_MM_ORDER[self.mm53] as usize;
                    let i5 = PARAM_FILE_MM_ORDER[self.mm35] as usize;
                    let i6 = PARAM_FILE_MM_ORDER[mm33 + 1] as usize;
                    tables.$field[i1][i2][i3][i4][i5][i6] = Some(val);
                }
                self.mm35 += 1;
                if self.mm35 == PARAM_FILE_MM_ORDER.len() {
                    self.mm35 = 1;
                    self.mm53 += 1;
                }
                if self.mm53 == PARAM_FILE_MM_ORDER.len() {
                    self.mm53 = 1;
                    self.mm55 += 1;
                }
                if self.mm55 == PARAM_FILE_MM_ORDER.len() {
                    self.mm55 = 1;
                    self.inner += 1;
                } 
                if self.inner == PARAM_FILE_PAIR_ORDER.len() - 1 {
                    self.inner = 0;
                    self.outer += 1;
                }
                if self.outer == PARAM_FILE_PAIR_ORDER.len() - 1 {
                    // The positions of pairtype: NN, and Base: N in the matrix.
                    let n_idx = [6usize, 6, 4, 4, 4, 4];
                    // Full ranges for each axis
                    let ranges: Vec<_> = n_idx.iter().map(|&n| 0..n).collect();
                    let masks: Vec<u8> = (1u8..64).collect();
                    for idx in ranges.iter().cloned().multi_cartesian_product() {
                        let val_idx: [usize; 6] = idx.clone().try_into().unwrap();
                        let val = tables.$field[val_idx[0]][val_idx[1]][val_idx[2]][val_idx[3]][val_idx[4]][val_idx[5]].unwrap();
                        //println!("{:?}: {}", val_idx, val);
                        for &mask in &masks {
                            //println!("{:06b}", mask);
                            let mut msk_idx = val_idx;
                            for axis in 0..6 {
                                if (mask >> axis) & 1 == 1 {
                                    msk_idx[axis] = n_idx[axis];
                                }
                            }
                            //println!("{:?}", val_idx);
                            //println!("{:?}", msk_idx);
                            tables.$field[msk_idx[0]][msk_idx[1]][msk_idx[2]][msk_idx[3]][msk_idx[4]][msk_idx[5]] = Some(
                                val.max(tables.$field[msk_idx[0]][msk_idx[1]][msk_idx[2]][msk_idx[3]][msk_idx[4]][msk_idx[5]].unwrap_or(i32::MIN)));
                        }
                    }
                }
            }
        }
    };
}

impl_int22_parser!(Int22, int22);
impl_int22_parser!(Int22Enthalpies, int22_enthalpies);

macro_rules! impl_loop_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name {
            base: usize,
        }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                let mut idx = 0;
                for token in line.split_whitespace()
                {
                    let val = if token == "INF" {
                        None
                    } else {
                        Some(token.parse::<i32>().expect(&format!(
                                    "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                                    stringify!($field), line, token)))
                    };
                    tables.$field[self.base+idx] = val;
                    idx += 1;
                }
                self.base += idx;
            }
        }
    };
}

impl_loop_parser!(Hairpin, hairpin);
impl_loop_parser!(HairpinEnthalpies, hairpin_enthalpies);
impl_loop_parser!(Bulge, bulge);
impl_loop_parser!(BulgeEnthalpies, bulge_enthalpies);
impl_loop_parser!(Interior, interior);
impl_loop_parser!(InteriorEnthalpies, interior_enthalpies);

macro_rules! impl_misc_parser {
    ($struct_name:ident, $field:ident, $ty:ty) => {
        #[derive(Default)]
        pub struct $struct_name;

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                let slice: Vec<i32> = line
                    .split_whitespace()
                    .map(|token| {
                        if let Ok(v) = token.parse::<i32>() {
                            v
                        } else if let Ok(f) = token.parse::<f64>() {
                            (f * 1000.0).round() as i32
                        } else {
                            panic!(
                                "Failed to parse number in {} while parsing line {:?}, token {:?}",
                                stringify!($field), line, token
                            )
                        }
                    })
                .collect();

                tables.$field = <$ty>::from_vrna_param_slice(&slice)
                    .expect(&format!("invalid values in {}", stringify!($field)));
            }
        }
    };
}

impl_misc_parser!(MLParams, ml_params, crate::nn_parsing::energy_tables::MLParams);
impl_misc_parser!(Ninio, ninio, crate::nn_parsing::energy_tables::NINIO);
impl_misc_parser!(Misc, misc, crate::nn_parsing::energy_tables::Misc);

macro_rules! impl_sequence_parser {
    ($struct_name:ident, $field:ident) => {
        #[derive(Default)]
        pub struct $struct_name { }

        impl SectionParser for $struct_name {
            fn parse_line(&mut self, line: &str, tables: &mut EnergyTables) {
                let mut parts = line.split_whitespace();
                if let (Some(seq), Some(g), Some(h)) = (parts.next(), parts.next(), parts.next()) {
                    let g = g.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            g
                    ));
                    let h = h.parse::<i32>().expect(&format!(
                            "Failed to parse integer in {} while parsing line {:?}, token {:?}",
                            stringify!($field),
                            line,
                            h
                    ));
                    tables.$field.insert(NucleotideVec::from_lossy(seq), (g, h));
                }
            }
        }
    };
}

impl_sequence_parser!(HairpinSequences, hairpin_sequences);


pub enum ParamFileSection {
    None,
    Stack(Stack),
    StackEnthalpies(StackEnthalpies),
    MismatchHairpin(MismatchHairpin),
    MismatchHairpinEnthalpies(MismatchHairpinEnthalpies),
    MismatchInterior(MismatchInterior),
    MismatchInteriorEnthalpies(MismatchInteriorEnthalpies),
    MismatchInterior1n(MismatchInterior1n),
    MismatchInterior1nEnthalpies(MismatchInterior1nEnthalpies),
    MismatchInterior23(MismatchInterior23),
    MismatchInterior23Enthalpies(MismatchInterior23Enthalpies),
    MismatchMulti(MismatchMulti),
    MismatchMultiEnthalpies(MismatchMultiEnthalpies),
    MismatchExterior(MismatchExterior),
    MismatchExteriorEnthalpies(MismatchExteriorEnthalpies),
    Dangle5(Dangle5),
    Dangle5Enthalpies(Dangle5Enthalpies),
    Dangle3(Dangle3),
    Dangle3Enthalpies(Dangle3Enthalpies),
    Int11(Int11),
    Int11Enthalpies(Int11Enthalpies),
    Int21(Int21),
    Int21Enthalpies(Int21Enthalpies),
    Int22(Int22),
    Int22Enthalpies(Int22Enthalpies),
    Hairpin(Hairpin),
    HairpinEnthalpies(HairpinEnthalpies),
    Bulge(Bulge),
    BulgeEnthalpies(BulgeEnthalpies),
    Interior(Interior),
    InteriorEnthalpies(InteriorEnthalpies),
    MLParams(MLParams),
    Ninio(Ninio),
    Misc(Misc),
    HairpinSequences(HairpinSequences),
}

macro_rules! section_match {
    ($s:expr, $($field:literal, $struct:ident),+ $(,)?) => {
        match $s {
            $(
                $field => Ok(ParamFileSection::$struct($struct::default())),
            )+
            _ => Err(()),
        }
    };
}

impl TryFrom<&str> for ParamFileSection {
    type Error = ();

    fn try_from(s: &str) -> Result<Self, ()> {
        let key = s.trim();
        section_match!(key,
            "stack", Stack, 
            "stack_enthalpies", StackEnthalpies, 
            "mismatch_hairpin", MismatchHairpin,
            "mismatch_hairpin_enthalpies", MismatchHairpinEnthalpies,
            "mismatch_interior", MismatchInterior,
            "mismatch_interior_enthalpies", MismatchInteriorEnthalpies,
            "mismatch_interior_1n", MismatchInterior1n,
            "mismatch_interior_1n_enthalpies", MismatchInterior1nEnthalpies,
            "mismatch_interior_23", MismatchInterior23,
            "mismatch_interior_23_enthalpies", MismatchInterior23Enthalpies,
            "mismatch_multi", MismatchMulti,
            "mismatch_multi_enthalpies", MismatchMultiEnthalpies,
            "mismatch_exterior", MismatchExterior,
            "mismatch_exterior_enthalpies", MismatchExteriorEnthalpies,
            "dangle5", Dangle5,
            "dangle5_enthalpies", Dangle5Enthalpies,
            "dangle3", Dangle3,
            "dangle3_enthalpies", Dangle3Enthalpies,
            "int11", Int11,
            "int11_enthalpies", Int11Enthalpies,
            "int21", Int21,
            "int21_enthalpies", Int21Enthalpies,
            "int22", Int22,
            "int22_enthalpies", Int22Enthalpies,
            "hairpin", Hairpin,
            "hairpin_enthalpies", HairpinEnthalpies,
            "bulge", Bulge,
            "bulge_enthalpies", BulgeEnthalpies,
            "interior", Interior,
            "interior_enthalpies", InteriorEnthalpies,
            "ML_params", MLParams,
            "NINIO", Ninio,
            "Misc", Misc,
            "Hexaloops", HairpinSequences,
            "Tetraloops", HairpinSequences,
            "Triloops", HairpinSequences,
        )
    }
}

