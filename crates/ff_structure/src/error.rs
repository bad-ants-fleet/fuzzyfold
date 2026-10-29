//! Errors for ff_structure. 
//!
//! NOTE: We communicate errors based on usize indexing, because error typically
//! occur when we have to cast between u16 <-> usize anyway.

use std::fmt;

/// Errors of parsing and converting structures. Positions are indices in
/// the structure string.
#[derive(Debug)]
pub enum StructureError {
    /// What was found, what was expected, and where.
    InvalidToken(String, String, usize),
    /// An opening bracket without a partner, at the first such position.
    UnmatchedOpen(char, usize),
    /// A closing bracket without a partner.
    UnmatchedClose(char, usize),
    /// The pair (i, j) needs more bracket levels than there are.
    TooManyPseudoknotLevels(usize, usize),
    /// The length of the structure, and the maximum length.
    TooLong(usize, usize),
}

impl fmt::Display for StructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            StructureError::UnmatchedOpen(c, i) | StructureError::UnmatchedClose(c, i) => {
                write!(f, "Unmatched '{}' at position {}", c, i)
            }
            StructureError::InvalidToken(tok, src, i) => {
                write!(f, "Invalid {} in {} at position {}", tok, src, i)
            }
            StructureError::TooManyPseudoknotLevels(i, j) => {
                write!(f, "Pair {} {} would require a larger pk-nesting support.", i, j)
            }
            StructureError::TooLong(n, max) => {
                write!(f, "Structure of length {} exceeds the maximum length {}", n, max)
            }
        }
    }
}

impl std::error::Error for StructureError {}


