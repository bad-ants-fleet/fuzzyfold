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
mod pair_list;

pub use representations::*;
pub use enumeration::*;
pub use dlfolding::*;
pub use pair_list::*;


