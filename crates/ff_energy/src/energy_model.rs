//! The `EnergyModel` trait and its errors.

use std::fmt;

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;

use crate::Pseudoknot;
use crate::RingDecomposition;
use crate::NestedLoop;
use crate::PairTypeRNA;
use crate::Base;

/// Errors of energy evaluation.
#[derive(Debug)]
pub enum EnergyError {
    /// A hairpin with fewer unpaired bases than the model allows.
    HairpinTooSmall { size: usize, min: usize },
    /// Two stacked pairs without parameters.
    UnsupportedStacking { outer: PairTypeRNA, inner: PairTypeRNA},
    /// A loop closed by bases that cannot pair.
    InvalidClosingPair,
    /// A circular strand together with other strands, in ViennaRNA.
    CircularStrand,
    /// A pseudoknot, in a model without pseudoknot parameters.
    UnsupportedPseudoknot,
    /// Sequence and structure differ in length.
    LengthMismatch { sequence: usize, structure: usize },
    /// Sequence and structure end a strand differently at this position.
    StrandMismatch(usize),
}

impl fmt::Display for EnergyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EnergyError::HairpinTooSmall { size, min } => {
                write!(f, "Hairpin too small: size {}, minimum allowed is {}", size, min)
            }
            EnergyError::UnsupportedStacking{ outer, inner } => {
                write!(f, "Stacking interaction [{}][{}] not supported by current energy model.", outer, inner)
            }
            EnergyError::InvalidClosingPair => {
                write!(f, "Invalid closing base pair")
            }
            EnergyError::CircularStrand => {
                write!(f, "ViennaRNA evaluates a circular strand only alone")
            }
            EnergyError::UnsupportedPseudoknot => {
                write!(f, "Pseudoknots are not supported by this energy model")
            }
            EnergyError::LengthMismatch { sequence, structure } => {
                write!(f, "Sequence length ({}) and structure length ({}) do not match", sequence, structure)
            }
            EnergyError::StrandMismatch(i) => {
                write!(f, "Sequence and structure end a strand differently at position {}", i)
            }
        }
    }
}

impl std::error::Error for EnergyError {}

/// A nearest neighbor energy model: the free energy of a structure is the
/// sum of the energies of its loops.
///
/// Energies are integers in units of 10 cal/mol (100 = 1 kcal/mol). The
/// sequence must have the length of the structure, with its strand breaks
/// and strand wraps at the same positions. This is not checked during
/// evaluation; see [`check_sequence`] for input. `Send + Sync`, so that simulations can share a model
/// between threads.
pub trait EnergyModel: Send + Sync {
    /// True if the two bases can form a pair.
    fn can_pair(&self, b1: Base, b2: Base) -> bool;

    /// The smallest number of unpaired bases in a hairpin.
    fn min_hairpin_size(&self) -> usize;

    /// The temperature in degrees Celsius.
    fn temperature(&self) -> f64;

    /// The free energy of a structure: the sum of its loops and
    /// pseudoknots. The first loop that cannot be evaluated ends the
    /// evaluation with its error.
    fn energy_of_structure<T: RingDecomposition>(&self,
        sequence: &[Base],
        structure: &T,
    ) -> Result<i32, EnergyError>;

    /// The free energy of a loop.
    fn energy_of_nested_loop(&self,
        sequence: &[Base],
        nested: &NestedLoop,
    ) -> Result<i32, EnergyError>;

    /// The free energy of a pseudoknot, with all its loops.
    /// `EnergyError::UnsupportedPseudoknot` in a model without pseudoknots.
    fn energy_of_pseudoknot(&self,
        sequence: &[Base],
        pseudoknot: &Pseudoknot,
    ) -> Result<i32, EnergyError>;

}

/// Checks that `structure` fits `sequence`: the same length, and strand
/// breaks and strand wraps at the same positions. Energy evaluation assumes
/// this without checking, so call this where both come from input.
pub fn check_sequence(sequence: &[Base], structure: &DotBracketVec) -> Result<(), EnergyError> {
    if sequence.len() != structure.len() {
        return Err(EnergyError::LengthMismatch {
            sequence: sequence.len(),
            structure: structure.len(),
        });
    }
    for (i, (&b, &d)) in sequence.iter().zip(structure.iter()).enumerate() {
        let fits = match b {
            Base::SB => d == DotBracket::STRAND_BREAK,
            Base::SW => d == DotBracket::STRAND_WRAP,
            _ => d != DotBracket::STRAND_BREAK && d != DotBracket::STRAND_WRAP,
        };
        if !fits {
            return Err(EnergyError::StrandMismatch(i));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::NucleotideVec;

    #[test]
    fn test_check_sequence() {
        let check = |seq: &str, db: &str| check_sequence(
            &NucleotideVec::try_from(seq).unwrap(),
            &DotBracketVec::try_from(db).unwrap(),
        );
        assert!(check("GGG+CCC", "(((+)))").is_ok());
        assert!(check("GGGAAACCC^", "(((...)))^").is_ok());
        assert!(matches!(check("GGGACCC", "(((.)))."), Err(EnergyError::LengthMismatch { .. })));
        assert!(matches!(check("GGG+CCC", "(((.)))"), Err(EnergyError::StrandMismatch(3))));
        assert!(matches!(check("GGG+CCC", "(((^)))"), Err(EnergyError::StrandMismatch(3))));
        assert!(matches!(check("GGGACCC", "(((+)))"), Err(EnergyError::StrandMismatch(3))));
    }
}
