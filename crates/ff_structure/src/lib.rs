//! Nucleic acid secondary structures: dot-bracket notation, pair tables,
//! and the closed regions of pseudoknotted structures.
//!
//! # Representations
//!
//! | Type               | Bytes per position | Used for                                        |
//! |--------------------|--------------------|-------------------------------------------------|
//! | [`DotBracketVec`]  | 1                  | storing, hashing and printing structures        |
//! | [`MultiPairTable`] | 2                  | energy evaluation, kinetics, structure analysis |
//!
//! Both hold one or more strands, with or without pseudoknots. Both parse
//! from a string with `TryFrom<&str>`, and convert into each other with
//! `TryFrom`. The end of a strand (`+`, `&` or `^`) takes one position, so
//! indices are positions in the structure string.
//!
//! ```
//! use ff_structure::{DotBracketVec, MultiPairTable, MultiStruct};
//!
//! let db = DotBracketVec::try_from("((.[[+))..]]").unwrap();
//! let mpt = MultiPairTable::try_from(&db).unwrap();
//! assert_eq!(mpt.get(0), MultiStruct::Paired(7));
//! assert_eq!(mpt.get(5), MultiStruct::StrandBreak);
//! assert_eq!(DotBracketVec::try_from(&mpt).unwrap(), db);
//! ```
//!
//! # Dot-bracket notation
//!
//! `.` is unpaired, `+` or `&` a strand break (written as `+`). Pairs use
//! `()`, `[]`, `{}`, `<>`, and the letters `Aa` to `Zz` (upper case opens):
//! 30 levels, see [`DotBracket`]. Each level matches on its own, so pairs of
//! different levels may cross.
//!
//! `^` is a strand wrap: the strand before it is circular, as in
//! `(((...)))^`. `ff_energy` evaluates circular strands, also together with
//! other strands and with pseudoknots; `ff_kinetics` simulates a circular
//! strand alone.
//!
//! Every strand has at least one position, so `+` and `^` must follow one,
//! and a structure does not end with `+`.
//!
//! Writing a structure from a pair table assigns the levels again, from
//! `(` upward, so a structure may come back with other brackets than it
//! was parsed with: `"([.)]"` stays as is, but `"[(.)]"` becomes `"((.))"`.
//!
//! # Pseudoknot analysis
//!
//! [`MultiPairTable::is_nested`] tells whether pairs cross. [`ClosedRegions`]
//! partitions a [`MultiPairTable`] into closed regions: nested regions,
//! closed by a pair, and pseudoknotted regions, made of bands; the density
//! counts the bands that span a common position.
//!
//! ```
//! use ff_structure::{ClosedRegions, MultiPairTable};
//!
//! let mpt = MultiPairTable::try_from("((..[[..))..]]").unwrap();
//! assert!(!mpt.is_nested());
//! let regions = ClosedRegions::from(&mpt);
//! assert_eq!(regions.iter().filter(|r| !r.nested).count(), 1);
//! assert_eq!(regions.density(), 2);
//! ```
//!
//! # Errors
//!
//! Parsing and conversions return a [`StructureError`], which carries the
//! position of the problem in the structure string.

mod error;
mod dotbracket;
mod multi_pair_table;
mod closed_regions;

pub use error::*;
pub use dotbracket::*;
pub use multi_pair_table::*;
pub use closed_regions::*;
