
use std::fmt;
use std::io::BufRead;
use ahash::AHashMap;

use crate::nn_parsing::ParamFileSection;
use crate::nn_parsing::SectionParser;
use crate::NucleotideVec;
use crate::BCOUNT as B;
use crate::PCOUNT as P;


#[derive(Debug)]
pub enum ParamError {
    Io(std::io::Error),
    Parse(String),
    MissingValue(&'static str, usize),
    InvalidLength(&'static str, usize, usize), 
    InvalidHairpinSize(usize),
}

impl std::error::Error for ParamError {}

impl From<std::io::Error> for ParamError {
    fn from(e: std::io::Error) -> Self {
        ParamError::Io(e)
    }
}

impl fmt::Display for ParamError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParamError::Io(e) => write!(f, "I/O error: {}", e),
            ParamError::Parse(msg) => write!(f, "Parse error: {}", msg),
            ParamError::MissingValue(table, index) => {
                write!(f, "Missing value in parameter table '{}' at index {}", table, index)
            }
            ParamError::InvalidLength(table, expected, got) => {
                write!(
                    f,
                    "Invalid length for parameter table '{}': expected {}, got {}",
                    table, expected, got
                )
            }
            ParamError::InvalidHairpinSize(n) => {
                write!(f, "Invalid hairpin size: {}", n)
            }
        }
    }
}

#[derive(Default)]
pub struct MLParams {
    pub base_en37: i32,
    pub base_enth: i32,
    pub closing_en37: i32,
    pub closing_enth: i32,
    pub intern_en37: i32,
    pub intern_enth: i32,
}

impl MLParams {
    pub fn from_vrna_param_slice(slice: &[i32]) -> Result<Self, ParamError> {
        if slice.len() != 6 {
            return Err(ParamError::InvalidLength("ml_params", 6, slice.len()));
        }
        Ok(Self {
            base_en37: slice[0],
            base_enth: slice[1],
            closing_en37: slice[2],
            closing_enth: slice[3],
            intern_en37: slice[4],
            intern_enth: slice[5],
        })
    }
}

#[derive(Default)]
pub struct NINIO {
    pub en37: i32,
    pub enth: i32,
    pub max: i32,
}

impl NINIO {
    pub fn from_vrna_param_slice(slice: &[i32]) -> Result<Self, ParamError> {
        if slice.len() != 3 {
            return Err(ParamError::InvalidLength("NINIO", 6, slice.len()));
        }
        Ok(Self {
            en37: slice[0],
            enth: slice[1],
            max: slice[2],
        })
    }
}

#[derive(Default)]
pub struct Misc {
    pub duplex_initiation_en37: i32,
    pub duplex_initiation_enth: i32,
    pub terminal_ru_en37: i32,
    pub terminal_ru_enth: i32,
    pub lxc: f64,
}

impl Misc {
    pub fn from_vrna_param_slice(slice: &[i32]) -> Result<Self, ParamError> {
        if slice.len() != 6 {
            // NOTE, we are discarding the last 0.
            return Err(ParamError::InvalidLength("Misc", 6, slice.len()));
        }
        Ok(Self {
            duplex_initiation_en37: slice[0],
            duplex_initiation_enth: slice[1],
            terminal_ru_en37: slice[2],
            terminal_ru_enth: slice[3],
            lxc: slice[4] as f64 / 1000.0,
        })
    }
}

fn rescale_energy(g_old: Option<i32>, h: Option<i32>, temp_change: f64) -> Option<i32> {
    match (g_old, h) {
        (Some(g), Some(h)) => {
            let g = (g) as f64;
            let h = (h) as f64;
            let s = h - g;
            Some((h - temp_change * s).round() as i32) 
            //Some((h - temp_change * s) as i32) // for better vrna compatibility, no rounding.
        }
        _ => None,
    }
}

trait RescaleWith {
    fn rescale_with(&mut self, enthalpies: &Self, temp_change: f64);
}

impl RescaleWith for Option<i32> {
    fn rescale_with(&mut self, enthalpies: &Self, temp_change: f64) {
        *self = rescale_energy(*self, *enthalpies, temp_change);
    }
}

impl<T: RescaleWith, const N: usize> RescaleWith for [T; N] {
    fn rescale_with(&mut self, enthalpies: &Self, temp_change: f64) {
        for (g, h) in self.iter_mut().zip(enthalpies.iter()) {
            g.rescale_with(h, temp_change);
        }
    }
}


pub struct EnergyTables {
    pub stack:            [[Option<i32>; P]; P],
    pub stack_enthalpies: [[Option<i32>; P]; P],
                                                                          
    pub mismatch_hairpin:            [[[Option<i32>; B]; B]; P],
    pub mismatch_hairpin_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_1n:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_1n_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_23:            [[[Option<i32>; B]; B]; P],
    pub mismatch_interior_23_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_multi:            [[[Option<i32>; B]; B]; P],
    pub mismatch_multi_enthalpies: [[[Option<i32>; B]; B]; P],
    pub mismatch_exterior:            [[[Option<i32>; B]; B]; P],
    pub mismatch_exterior_enthalpies: [[[Option<i32>; B]; B]; P],
                                                                          
    pub dangle5:            [[Option<i32>; B]; P],
    pub dangle5_enthalpies: [[Option<i32>; B]; P],
    pub dangle3:            [[Option<i32>; B]; P],
    pub dangle3_enthalpies: [[Option<i32>; B]; P],
                                                                          
    pub int11:            Box<[[[[Option<i32>; B]; B]; P]; P]>,
    pub int11_enthalpies: Box<[[[[Option<i32>; B]; B]; P]; P]>,
    pub int21:            Box<[[[[[Option<i32>; B]; B]; B]; P]; P]>,
    pub int21_enthalpies: Box<[[[[[Option<i32>; B]; B]; B]; P]; P]>,
    pub int22:            Box<[[[[[[Option<i32>; B]; B]; B]; B]; P]; P]>,
    pub int22_enthalpies: Box<[[[[[[Option<i32>; B]; B]; B]; B]; P]; P]>,

    pub hairpin:            [Option<i32>; 31],
    pub hairpin_enthalpies: [Option<i32>; 31],
    pub bulge:            [Option<i32>; 31],
    pub bulge_enthalpies: [Option<i32>; 31],
    pub interior:            [Option<i32>; 31],
    pub interior_enthalpies: [Option<i32>; 31],

    pub ml_params: MLParams,
    pub ninio: NINIO,
    pub misc: Misc,

    pub hairpin_sequences: AHashMap<NucleotideVec, (i32, i32)>,
}

impl Default for EnergyTables {
    fn default() -> Self {
        Self {
            stack:            [[None; P]; P],
            stack_enthalpies: [[None; P]; P],

            mismatch_hairpin:            [[[None; B]; B]; P],
            mismatch_hairpin_enthalpies: [[[None; B]; B]; P],
            mismatch_interior:            [[[None; B]; B]; P],
            mismatch_interior_enthalpies: [[[None; B]; B]; P],
            mismatch_interior_1n:            [[[None; B]; B]; P],
            mismatch_interior_1n_enthalpies: [[[None; B]; B]; P],
            mismatch_interior_23:            [[[None; B]; B]; P],
            mismatch_interior_23_enthalpies: [[[None; B]; B]; P],
            mismatch_multi:            [[[None; B]; B]; P],
            mismatch_multi_enthalpies: [[[None; B]; B]; P],
            mismatch_exterior:            [[[None; B]; B]; P],
            mismatch_exterior_enthalpies: [[[None; B]; B]; P],
            dangle5:            [[None; B]; P],
            dangle5_enthalpies: [[None; B]; P],
            dangle3:            [[None; B]; P],
            dangle3_enthalpies: [[None; B]; P],

            int11:            Box::new([[[[None; B]; B]; P]; P]),
            int11_enthalpies: Box::new([[[[None; B]; B]; P]; P]),
            int21:            Box::new([[[[[None; B]; B]; B]; P]; P]),
            int21_enthalpies: Box::new([[[[[None; B]; B]; B]; P]; P]),
            int22:            Box::new([[[[[[None; B]; B]; B]; B]; P]; P]),
            int22_enthalpies: Box::new([[[[[[None; B]; B]; B]; B]; P]; P]),

            hairpin: [None; 31],
            hairpin_enthalpies: [None; 31],
            bulge: [None; 31],
            bulge_enthalpies: [None; 31],
            interior: [None; 31],
            interior_enthalpies: [None; 31],
            ml_params: MLParams::default(),
            ninio: NINIO::default(),
            misc: Misc::default(),

            hairpin_sequences: AHashMap::new(),
        }
    }
}

macro_rules! section_match {
    ($enum:expr, $line:expr, $tables:expr, $($struct:ident),+ $(,)?) => {
        match $enum {
            $(
                ParamFileSection::$struct(ref mut s) => s.parse_line($line, &mut $tables),
            )+
                _ => { panic!("Couldn't match line \"{}\"", $line) },
        }
    };
}

impl EnergyTables {
    pub fn from_reader<R: BufRead>(reader: R) -> Result<Self, ParamError> {
        let mut tables = EnergyTables::default();
        let mut section = ParamFileSection::None;

        for line in reader.lines() {
            let line = line?;
            let line = line.trim();

            if line.is_empty() || line.starts_with("/*"){
                continue;
            } 

            if let Some(rest) = line.strip_prefix("# ") {
                match ParamFileSection::try_from(rest.trim()) {
                    Ok(sec) => {
                        section = sec;
                    }
                    Err(_) => {
                        eprintln!("⚠️  Unknown parameter file section: {:?}", rest);
                        return Err(ParamError::Parse(rest.to_string()));
                    }
                }
                continue;
            } else if line.starts_with("#") {
                continue;
            }

            section_match!(section, line, tables, 
                Stack, 
                StackEnthalpies,
                MismatchHairpin,
                MismatchHairpinEnthalpies,
                MismatchInterior,
                MismatchInteriorEnthalpies,
                MismatchInterior1n,
                MismatchInterior1nEnthalpies,
                MismatchInterior23,
                MismatchInterior23Enthalpies,
                MismatchMulti,
                MismatchMultiEnthalpies,
                MismatchExterior,
                MismatchExteriorEnthalpies,
                Dangle5,
                Dangle5Enthalpies,
                Dangle3,
                Dangle3Enthalpies,
                Int11,
                Int11Enthalpies,
                Int21,
                Int21Enthalpies,
                Int22,
                Int22Enthalpies,
                Hairpin,
                HairpinEnthalpies,
                Bulge,
                BulgeEnthalpies,
                Interior,
                InteriorEnthalpies,
                MLParams,
                Ninio,
                Misc,
                HairpinSequences,
            );
        }

        Ok(tables) 
    }

    pub fn rescale(&mut self, temp_change: f64) {
        self.stack.rescale_with(&self.stack_enthalpies, temp_change);
        self.mismatch_hairpin.rescale_with(&self.mismatch_hairpin_enthalpies, temp_change);
        self.mismatch_interior.rescale_with(&self.mismatch_interior_enthalpies, temp_change);
        self.mismatch_interior_1n.rescale_with(&self.mismatch_interior_1n_enthalpies, temp_change);
        self.mismatch_interior_23.rescale_with(&self.mismatch_interior_23_enthalpies, temp_change);
        self.mismatch_multi.rescale_with(&self.mismatch_multi_enthalpies, temp_change);
        self.mismatch_exterior.rescale_with(&self.mismatch_exterior_enthalpies, temp_change);
        self.dangle5.rescale_with(&self.dangle5_enthalpies, temp_change);
        self.dangle3.rescale_with(&self.dangle3_enthalpies, temp_change);
        (*self.int11).rescale_with(&*self.int11_enthalpies, temp_change);
        (*self.int21).rescale_with(&*self.int21_enthalpies, temp_change);
        (*self.int22).rescale_with(&*self.int22_enthalpies, temp_change);
        self.hairpin.rescale_with(&self.hairpin_enthalpies, temp_change);
        self.bulge.rescale_with(&self.bulge_enthalpies, temp_change);
        self.interior.rescale_with(&self.interior_enthalpies, temp_change);

        self.ml_params.base_en37 = rescale_energy(
            Some(self.ml_params.base_en37),
            Some(self.ml_params.base_enth),
            temp_change).unwrap();

        self.ml_params.closing_en37 = rescale_energy(
            Some(self.ml_params.closing_en37),
            Some(self.ml_params.closing_enth),
            temp_change).unwrap();

        self.ml_params.intern_en37 = rescale_energy(
            Some(self.ml_params.intern_en37),
            Some(self.ml_params.intern_enth),
            temp_change).unwrap();

        self.ninio.en37 = rescale_energy(
            Some(self.ninio.en37),
            Some(self.ninio.enth),
            temp_change).unwrap();

        self.misc.duplex_initiation_en37 = rescale_energy(
            Some(self.misc.duplex_initiation_en37),
            Some(self.misc.duplex_initiation_enth),
            temp_change).unwrap();

        self.misc.terminal_ru_en37 = rescale_energy(
            Some(self.misc.terminal_ru_en37),
            Some(self.misc.terminal_ru_enth),
            temp_change).unwrap();

        self.misc.lxc *= temp_change;

        self.hairpin_sequences.iter_mut()
            .for_each(|(_key, (g, h))| {
                *g = rescale_energy(Some(*g), Some(*h), temp_change).unwrap();
            });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use crate::PairTypeRNA;
    use crate::Base;

    #[test]
    fn test_parse_stack() {
        let dummy = r#"
# stack
/*  CG    GC    GU    UG    AU    UA    NN          */
  -240  -330  -210  -140  -210  -210  -140    /* CG */
  -330  -340  -250  -150  -220  -240  -150    /* GC */
  -210  -250   130   -50  -140  -130   130    /* GU */
  -140  -150   -50    30   -60  -100    30    /* UG */
  -210  -220  -140   -60  -110   -90   -60    /* AU */
  -210  -240  -130  -100   -90  -130   -90    /* UA */
  -140  -150   130    30   -60   -90   130    /* NN */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.stack[PairTypeRNA::CG as usize][PairTypeRNA::CG as usize], Some(-240));
        assert_eq!(tables.stack[PairTypeRNA::GC as usize][PairTypeRNA::CG as usize], Some(-330));
        assert_eq!(tables.stack[PairTypeRNA::GU as usize][PairTypeRNA::CG as usize], Some(-210));
    }

    #[test]
    fn test_parse_mismatch() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# mismatch_hairpin                                                                                         
  -80  -100  -110  -100   -80    /* CG,E */   
 -140  -150  -150  -140  -150    /* CG,A */
  -80  -100  -110  -100   -80    /* CG,C */   
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.mismatch_hairpin[CG as usize][N as usize][N as usize], Some(-80));
        assert_eq!(tables.mismatch_hairpin[CG as usize][N as usize][A as usize], Some(-100));
    }

    #[test]
    fn test_parse_dangle() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# dangle5
/*   N     A     C     G     U          */
   -10   -50   -30   -20   -10    /* CG */
    -0   -20   -30    -0    -0    /* GC */
   -20   -30   -30   -40   -20    /* GU */
   -10   -30   -10   -20   -20    /* UG */
   -20   -30   -30   -40   -20    /* AU */
   -10   -30   -10   -20   -20    /* UA */
    -0   -20   -10    -0    -0    /* NN */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.dangle5[CG as usize][N as usize], Some(-10));
        assert_eq!(tables.dangle5[CG as usize][A as usize], Some(-50));
    }

    #[test]
    fn test_parse_int11() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# int11
  90    90    50    50    50    /* CG,CG,N */
  90    90    50    50    50    /* CG,CG,A */
  50    50    50    50    50    /* CG,CG,C */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][N as usize], Some(90));
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][A as usize], Some(90));
        assert_eq!(tables.int11[CG as usize][CG as usize][N as usize][C as usize], Some(50));
    }

    #[test]
    fn test_parse_int21() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"
# int21
   230   230   230   230   230    /* CG,CG,N,N */
   230   230   230   230   230    /* CG,CG,N,A */
   230   230   230   230   230    /* CG,CG,N,C */
   230   230   230   230   230    /* CG,CG,N,G */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.int21[CG as usize][CG as usize][N as usize][N as usize][N as usize], Some(230));
        assert_eq!(tables.int21[CG as usize][CG as usize][N as usize][N as usize][A as usize], Some(230));
    }

    #[test]
    fn test_parse_int22() {
        use Base::*;
        use PairTypeRNA::*;
        let dummy = r#"

# int22
   120   160    20   160    /* CG,CG,A,A,A */
   110   150    20   150    /* CG,CG,A,A,C */
    20    60   -70    60    /* CG,CG,A,A,G */
   110   150    20   150    /* CG,CG,A,A,U */
   160   200    60   200    /* CG,CG,A,C,A */
   140   180   110   180    /* CG,CG,A,C,C */
   160   200    60   200    /* CG,CG,A,C,G */
   130   170    90   170    /* CG,CG,A,C,U */
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.int22[CG as usize][CG as usize][A as usize][A as usize][A as usize][A as usize], Some(120));
        assert_eq!(tables.int22[CG as usize][CG as usize][A as usize][A as usize][A as usize][C as usize], Some(160));
    }


    #[test]
    fn test_parse_loops() {
        let dummy = r#"
# hairpin
   INF   INF   INF   540   560   570   540   600   550   640
   650   660   670   680   690   690   700   710   710   720
   720   730   730   740   740   750   750   750   760   760
   770
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();

        // Check one of the parsed entries
        assert_eq!(tables.hairpin[0], None);
        assert_eq!(tables.hairpin[1], None);
        assert_eq!(tables.hairpin[2], None);
        assert_eq!(tables.hairpin[3], Some(540));
        assert_eq!(tables.hairpin[29], Some(760));
        assert_eq!(tables.hairpin[30], Some(770));
    }

    #[test]
    fn test_sequence_parsing() {
        let dummy = r#"
# Hexaloops
ACAGUACU     280   -1680
ACAGUGAU     360   -1140
ACAGUGCU     290   -1280
ACAGUGUU     180   -1540

# Tetraloops
CAACGG     550     690
CCAAGG     330   -1030
CCACGG     370    -330
CCCAGG     340    -890
CCGAGG     350    -660
CCGCGG     360    -750
CCUAGG     370    -350
CCUCGG     250   -1390
CUAAGG     360    -760
CUACGG     280   -1070
CUCAGG     370    -660
CUCCGG     270   -1290
CUGCGG     280   -1070
CUUAGG     350    -620
CUUCGG     370   -1530
CUUUGG     370    -680

# Triloops
CAACG     680    2370
GUUAC     690    1080
"#;

        let cursor = Cursor::new(dummy);
        let tables = EnergyTables::from_reader(cursor).unwrap();
        println!("{:?}", tables.hairpin_sequences);

        // Check one of the parsed entries
        assert_eq!(tables.hairpin_sequences[&NucleotideVec::from_lossy("CCAAGG")], (330, -1030));
        assert_eq!(tables.hairpin_sequences[&NucleotideVec::from_lossy("CAACG")], (680, 2370));
    }

}

