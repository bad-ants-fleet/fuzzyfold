//! Errors for ff_structure. 
//!
//! NOTE: We communicate errors based on usize indexing, because error typically
//! occur when we have to cast between u16 <-> usize anyway.

use std::fmt;


/// Errors that can occur when parsing a structure string into a PairTable
/// different variants with associated data for different types of errors:
#[derive(Debug)]
pub enum StructureError {
    InvalidToken(String, String, usize),
    UnmatchedOpen(usize),
    UnmatchedClose(usize),
    UnmatchedMultiOpen((usize, usize)),
    UnmatchedMultiClose((usize, usize)),
}

impl fmt::Display for StructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureError::UnmatchedOpen(i) => {
                write!(f, "Unmatched '(' at position {}", i)
            }
            StructureError::UnmatchedMultiOpen((si, di)) => {
                write!(f, "Unmatched '(' at strand {}, domain {}", si, di)
            }
            StructureError::UnmatchedClose(i) => {
                write!(f, "Unmatched ')' at position {}", i)
            }
            StructureError::UnmatchedMultiClose((si, di)) => {
                write!(f, "Unmatched ')' at strand {}, domain {}", si, di)
            }
            StructureError::InvalidToken(tok, src, i) => {
                write!(f, "Invalid {} in {} at position {}", tok, src, i)
            }
        }
    }
}

impl std::error::Error for StructureError {}


