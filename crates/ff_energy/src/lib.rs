//! Free energies of nucleic acid secondary structures in the nearest neighbor
//! model, including pseudoknots and structures of several strands.
//!
//! A structure is split into loops traced along the backbone and
//! pseudoknots ([`RingDecomposition`], implemented for
//! `ff_structure::MultiPairTable`), and an [`EnergyModel`] sums their
//! energies. A [`Pseudoknot`] is not a loop: it owns the loops on cycles of
//! the loop graph (a loop per vertex, a pair per edge), and which of their
//! segments are strained. Two models implement it:
//! - [`Fuzzyfold`]: the model of fuzzyfold, under development. The loops of
//!   ViennaRNA, with pseudoknot terms on top ([`parameters::PKParams`]);
//!   the open circle costs the initiation of a hairpin loop of its ring.
//! - [`ViennaRNA`]: the model of ViennaRNA 2.6, without pseudoknots, kept
//!   as it is for comparisons.
//!
//! ```
//! use ff_structure::MultiPairTable;
//! use ff_energy::{EnergyModel, NucleotideVec, Fuzzyfold};
//! use ff_energy::parameters::RNA_TURNER_2004;
//!
//! let model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
//! let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
//! let pk = MultiPairTable::try_from("[[[....(((((]]]......)))))").unwrap();
//! assert_eq!(model.energy_of_structure(&seq, &pk).unwrap(), -357);
//! ```
//!
//! # Units
//!
//! Energies are integers in units of 10 cal/mol: -357 is -3.57 kcal/mol.
//! Temperatures are in degrees Celsius.
//!
//! # Sequences
//!
//! A [`NucleotideVec`] holds the bases A, C, G, U, T, pseudouridine (`P`) and
//! inosine (`I`), strand breaks (`+` or `&`) and strand wraps (`^`, the strand
//! before it is circular). The end of a strand takes one position, as in the
//! structure, so that indices agree.
//!
//! # Parameter sets
//!
//! In [`parameters`]:
//! - `RNA_TURNER_2004`: RNA, Turner 2004.
//! - `RNA_EXTENDED`: Turner 2004 with special hairpins and modified bases.
//!   The default of `Fuzzyfold::default()`.
//! - `RNA_ANDRONESCU_2007`: RNA, fitted parameters, at 37 C only
//!   (`Fuzzyfold::from_andronescu_params`).
//! - `DNA_MATHEWS_2004`: DNA.
//!
//! `Fuzzyfold::from_thermo_params(&set, celsius)` rescales a set to the
//! temperature. The pseudoknot parameters are the same for every set, given
//! at 37 C and not rescaled; using them at another temperature logs a warning.
//!
//! # Loops and terms
//!
//! [`Fuzzyfold::energy_terms`] lists the energy of every loop, and for a
//! pseudoknot its terms, as `ff-eval -v` prints them. The energies of
//! level 0 sum to the energy of the structure.
//!
//! ```
//! use ff_structure::MultiPairTable;
//! use ff_energy::{Ring, RingDecomposition, NucleotideVec, Fuzzyfold};
//! use ff_energy::parameters::RNA_TURNER_2004;
//!
//! let model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
//! let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
//! let pk = MultiPairTable::try_from("[[[....(((((]]]......)))))").unwrap();
//!
//! // The stacks of the stems, and the pseudoknot: the exterior loop is
//! // part of its pseudoloop.
//! let rings = pk.rings();
//! assert!(matches!(rings.last(), Some(Ring::Pseudoknot(_))));
//!
//! let terms = model.energy_terms(&seq, &pk).unwrap();
//! let total: i32 = terms.iter().filter(|(level, _, _)| *level == 0).map(|t| t.2).sum();
//! assert_eq!(total, -357);
//! ```
//!
//! # Errors
//!
//! Parsing a sequence returns a [`SequenceError`], evaluating a structure an
//! [`EnergyError`].

/// Parameter sets.
pub mod parameters;

mod models;
mod nucleotides;
mod decomposition;
mod energy_model;

pub use nucleotides::*;
pub use decomposition::*;
pub use energy_model::*;
pub use models::*;
