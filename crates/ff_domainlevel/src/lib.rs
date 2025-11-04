//! The ff_domainlevel crate.
//!
//! Provides basic domain-level representations of nucleic acids:
//!  - Domains
//!  - Strands
//!  - Complexes
//!  - Reactions
//!
//! Provides some domain-level folding utilities.
//!  - base-pair maximization (Nussinov).
//!  - domain-level reaction enumeration.
//!  - PairList representatin (1-based).
//!

/// Design module, mostly ACFP stuff.
pub mod design;

/// Rules module, used by enumeration.
pub mod rules;

mod representations;
mod enumeration;
mod dlfolding;
mod pair_set;

pub use representations::*;
pub use enumeration::*;
pub use dlfolding::*;
pub use pair_set::*;

pub type P1KEY = u32;


