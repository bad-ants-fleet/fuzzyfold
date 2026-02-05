use std::fs::File;
use std::io::{stdin, BufRead, BufReader, Cursor};
use std::path::Path;

use anyhow::{anyhow, Result};
use paste::paste;
use ff_structure::DotBracketVec;
use ff_energy::NucleotideVec;

// ============================================================
//  Generic FASTA-like parser supporting lenient/strict modes
// ============================================================

#[derive(Clone, Copy)]
enum FastaMode {
    Lenient,
    Strict,
}

/// Core parsing logic shared by all adapters.
fn parse_fasta_like<R: BufRead>(
    reader: R,
    mode: FastaMode,
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
            sequence = Some(NucleotideVec::from_lossy(token));
            //sequence = Some(NucleotideVec::try_from(token)?);
        } else if structure.is_none() {
            let token = line.split_whitespace().next().unwrap();
            structure = Some(DotBracketVec::try_from(token)?);
            break;
        }
    }

    let sequence = sequence.ok_or_else(|| anyhow!("Missing sequence line"))?;

    let structure = match (structure, mode) {
        (Some(s), FastaMode::Strict) => {
            if sequence.len() != s.len() {
                return Err(anyhow!(
                        "Sequence length ({}) and structure length ({}) do not match",
                        sequence.len(),
                        s.len()
                ));
            }
            s
        },
        (None, FastaMode::Strict) => return Err(anyhow!("Missing structure line")),

        (Some(s), FastaMode::Lenient) => {
            if sequence.len() < s.len() {
                return Err(anyhow!(
                        "Structure is longer than sequence ({} > {}).",
                        s.len(), sequence.len()
                ));
            }
            s
        },

        (None, FastaMode::Lenient) => {
            DotBracketVec::try_from(".")
                .expect("Failed to construct open-chain structure")
        }
    };

    Ok((header, sequence, structure))
}

// ============================================================
//  Base parser functions (lenient and strict variants)
// ============================================================

pub fn read_cotr<R: BufRead>(reader: R) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    parse_fasta_like(reader, FastaMode::Lenient)
}

pub fn read_eval<R: BufRead>(reader: R) -> Result<(Option<String>, NucleotideVec, DotBracketVec)> {
    parse_fasta_like(reader, FastaMode::Strict)
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
/// define_input_variants!(read_fasta_like, Result<(Option<String>, NucleotideVec, DotBracketVec)>);
/// ```
macro_rules! define_input_variants {
    ($base:ident, $ret:ty) => {
        paste! {
            /// Read from a string buffer.
            pub fn [<$base _string>](s: &str) -> $ret {
                $base(Cursor::new(s))
            }

            /// Read from a file path.
            pub fn [<$base _file>]<P: AsRef<Path>>(path: P) -> $ret {
                let reader = BufReader::new(File::open(path)?);
                $base(reader)
            }

            /// Read from stdin.
            pub fn [<$base _stdin>]() -> $ret {
                let reader = BufReader::new(stdin());
                $base(reader)
            }

            /// Read either from stdin ("-") or a file path.
            pub fn [<$base _input>](s: &str) -> $ret {
                if s == "-" {
                    [<$base _stdin>]()
                } else {
                    [<$base _file>](s)
                }
            }
        }
    };
}

// ============================================================
//  Apply macro to generate adapters for both variants
// ============================================================

type FastaResult = Result<(Option<String>, NucleotideVec, DotBracketVec)>;

define_input_variants!(read_cotr, FastaResult);
define_input_variants!(read_eval, FastaResult);

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
        let (hdr, seq, dbv) = read_cotr_string(input).unwrap();
        assert_eq!(hdr, Some(">test".into()));
        assert_eq!(seq.to_string(), "ACGU");
        assert_eq!(dbv.to_string(), "....");

        let input = ">test\nACGU";
        let (hdr, seq, dbv) = read_cotr_string(input).unwrap();
        assert_eq!(hdr, Some(">test".into()));
        assert_eq!(seq.to_string(), "ACGU");
        assert_eq!(dbv.to_string(), ".");
    }

    #[test]
    fn test_read_eval_input() {
        let input = ">test\nACGU\n....\n";
        let ok = read_eval_string(input);
        assert!(ok.is_ok());

        let missing = ">test\nACGU\n";
        let err = read_eval_string(missing);
        assert!(err.is_err(), "Missing structure line should fail in strict mode");
    }
}

