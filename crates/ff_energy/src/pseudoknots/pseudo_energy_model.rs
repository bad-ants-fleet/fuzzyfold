//! [`PseudoEnergyModel`] trait for evaluating energies of pseudoknotted structures.

use super::{ClosingDescriptor, Loop};
use crate::{Base, EnergyError};

/// Energy model for pseudoknotted RNA structures.
///
/// Evaluates a [`Vec<Loop>`] produced by [`crate::parse_structure`].
/// Standard loop types (Stack, Hairpin, Interior, Bulge, Multiloop, External)
/// delegate to an underlying nearest-neighbor model.  [`LoopType::Pseudoloop`]
/// uses a simplified multiloop approximation:
///
/// ```text
/// E_pseudoloop ≈ ml_closing + ml_intern × (n_children + 2)
/// ```
///
/// where the two closing pairs of the pseudoloop each count as one stem.
/// Unpaired-base and mismatch contributions are omitted for the pseudoloop —
/// those are deferred to a future parameter extension.
pub trait PseudoEnergyModel {
    fn energy_of_pseudo_loop(&self, sequence: &[Base], lp: &Loop) -> Result<i32, EnergyError>;

    fn energy_of_pseudoknotted_structure(
        &self,
        sequence: &[Base],
        loops: &[Loop],
    ) -> Result<i32, EnergyError> {
        loops.iter().try_fold(0i32, |acc, lp| {
            self.energy_of_pseudo_loop(sequence, lp).map(|e| acc + e)
        })
    }
}

/// Helper: extract the single pair from a `ClosingDescriptor::Single`, or error.
pub fn single_pair(cd: Option<ClosingDescriptor>) -> Result<(usize, usize), EnergyError> {
    match cd {
        Some(ClosingDescriptor::Single(p)) => Ok(p),
        _ => Err(EnergyError::InvalidClosingPair),
    }
}

type TwoClosingPairs = ((usize, usize), (usize, usize));

/// Helper: extract both pairs from a `ClosingDescriptor::Double`, or error.
pub fn double_pair(closing_desc: Option<ClosingDescriptor>) -> Result<TwoClosingPairs, EnergyError> {
    match closing_desc {
        Some(ClosingDescriptor::Double(p1, p2)) => Ok((p1, p2)),
        _ => Err(EnergyError::InvalidClosingPair),
    }
}

/// Helper: collect branch pairs from children, failing on any `Double` descriptor.
pub fn collect_single_branches(
    children: &[ClosingDescriptor],
) -> Result<Vec<(usize, usize)>, EnergyError> {
    children
        .iter()
        .map(|cd| match cd {
            ClosingDescriptor::Single(p) => Ok(*p),
            ClosingDescriptor::Double(..) => Err(EnergyError::InvalidClosingPair),
        })
        .collect()
}
