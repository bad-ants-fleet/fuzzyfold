use std::fs::File;
use std::io::{stdin, BufRead, BufReader, Cursor};
use std::path::Path;

use anyhow::{anyhow, Result};
use paste::paste;
use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;

#[derive(Clone, Copy)]
enum NAMode {
    Lenient,
    Strict,
}

/// Core parsing logic shared by all adapters.
fn parse_na_format<R: BufRead>(
    reader: R,
    mode: NAMode,
    is_rna: bool,
) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    let mut header: Option<String> = None;
    let mut sequence: Option<NucleotideVec> = None;
    let mut structure: Option<DotBracketVec> = None;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            if sequence.is_some() && structure.is_some() {
                break;
            } else {
                continue;
            }
        }

        if line.starts_with('>') {
            header = Some(line.to_string());
        } else if sequence.is_none() {
            let token = line.split_whitespace().next().unwrap();
            if is_rna {
                sequence = Some(NucleotideVec::try_from_rna(token)?);
            } else {
                sequence = Some(NucleotideVec::try_from_dna(token)?);
            }
        } else if structure.is_none() {
            let token = line.split_whitespace().next().unwrap();
            structure = Some(DotBracketVec::try_from(token)?);
            break;
        }
    }

    let sequence = sequence.ok_or_else(|| anyhow!("Missing sequence line"))?;

    let structure = match (structure, mode) {
        (Some(s), NAMode::Strict) => {
            if sequence.len() != s.len() {
                return Err(anyhow!(
                        "Sequence length ({}) and structure length ({}) do not match",
                        sequence.len(),
                        s.len()
                ));
            }
            s
        },
        (None, NAMode::Strict) => return Err(anyhow!("Missing structure line")),

        (Some(s), NAMode::Lenient) => {
            if sequence.len() < s.len() {
                return Err(anyhow!(
                        "Structure is longer than sequence ({} > {}).",
                        s.len(), sequence.len()
                ));
            }
            s
        },

        (None, NAMode::Lenient) => {
            DotBracketVec::try_from(".")
                .expect("Failed to construct open-chain structure")
        }
    };

    Ok((header, sequence, structure))
}

fn read_fasta<R: BufRead>(reader: R, _rna: bool
) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    let mut header: Option<String> = None;
    let mut seq_chars = String::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            if !seq_chars.is_empty() { break; }
            continue;
        }
        if line.starts_with('>') {
            if seq_chars.is_empty() {
                header = Some(line.to_string());
            } else {
                break; // next record
            }
        } else {
            seq_chars.push_str(line);
        }
    }

    if seq_chars.is_empty() {
        return Err(anyhow!("Missing sequence"));
    }

    let sequence = NucleotideVec::try_from_dna(&seq_chars)?;
    let structure = DotBracketVec::try_from(".").expect("open-chain fallback");

    Ok((header, sequence, structure))
}

// ============================================================
//  Base parser functions (lenient and strict variants)
// ============================================================

pub fn read_cotr<R: BufRead>(reader: R, is_rna: bool) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    parse_na_format(reader, NAMode::Lenient, is_rna)
}

pub fn read_eval<R: BufRead>(reader: R, is_rna: bool) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    parse_na_format(reader, NAMode::Strict, is_rna)
}

// ============================================================
//  Macro generating file/string/stdin/input helpers
// ============================================================

/// Generate input adapters for a base parser function `fn base<R: BufRead>(R) -> Result<T>`.
///
/// This expands into:
/// - `base_string(&str)`
/// - `base_file<P: AsRef<Path>>(P)`
/// - `base_stdin()`
/// - `base_input(&str)`  (dispatches "-" → stdin, otherwise → file)
///
/// Example:
/// ```ignore
/// define_input_variants!(read_na_format, Result<(Option<String>, NucleotideVec, DotBracketVec)>);
/// ```
macro_rules! define_input_variants {
    ($base:ident, $ret:ty) => {
        paste! {
            /// Read from a string buffer.
            pub fn [<$base _string>](s: &str, rna: bool) -> $ret {
                $base(Cursor::new(s), rna)
            }

            /// Read from a file path.
            pub fn [<$base _file>]<P: AsRef<Path>>(path: P, rna: bool) -> $ret {
                let reader = BufReader::new(File::open(path)?);
                $base(reader, rna)
            }

            /// Read from stdin.
            pub fn [<$base _stdin>](rna: bool) -> $ret {
                let reader = BufReader::new(stdin());
                $base(reader, rna)
            }

            /// Read either from stdin ("-") or a file path.
            pub fn [<$base _input>](s: &str, rna: bool) -> $ret {
                if s == "-" {
                    [<$base _stdin>](rna)
                } else {
                    [<$base _file>](s, rna)
                }
            }
        }
    };
}

// ============================================================
//  Apply macro to generate adapters for both variants
// ============================================================

type NAResult = Result<(Option<String>, NucleotideVec, DotBracketVec)>;

define_input_variants!(read_cotr, NAResult);
define_input_variants!(read_eval, NAResult);
define_input_variants!(read_fasta, NAResult);

// ============================================================
//  Example helper: ruler()
// ============================================================

pub fn ruler(len: usize) -> String {
    let mut s = String::new();
    let mut c = 0;
    for i in 0..=len {
        if i % 10 == 0 {
            let t = format!("{}", i / 10);
            c = t.len() - 1;
            s.push_str(&t);
            continue;
        } else if c > 0 {
            c -= 1;
            continue;
        }
        if i % 10 == 5 {
            s.push(',');
        } else {
            s.push('.');
        }
    }
    s
}

// ============================================================
//  Unit tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ruler() {
        assert_eq!(ruler(0), "0");
        assert_eq!(ruler(5), "0....,");
        assert_eq!(ruler(10), "0....,....1");
    }

    #[test]
    fn test_read_cotr_input() {
        let input = ">test\nACGU\n....\n";
        let (hdr, seq, dbv) = read_cotr_string(input, true).unwrap();
        assert_eq!(hdr, Some(">test".into()));
        assert_eq!(seq.to_string(), "ACGU");
        assert_eq!(dbv.to_string(), "....");

        let input = ">test\nACGU";
        let (hdr, seq, dbv) = read_cotr_string(input, true).unwrap();
        assert_eq!(hdr, Some(">test".into()));
        assert_eq!(seq.to_string(), "ACGU");
        assert_eq!(dbv.to_string(), ".");
    }

    #[test]
    fn test_read_eval_input() {
        let input = ">test\nACGU\n....\n";
        let ok = read_eval_string(input, true);
        assert!(ok.is_ok());

        let missing = ">test\nACGU\n";
        let err = read_eval_string(missing, true);
        assert!(err.is_err(), "Missing structure line should fail in strict mode");
    }
}

