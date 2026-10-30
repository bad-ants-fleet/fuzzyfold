
use ff_structure::DotBracketVec;
use ff_energy::EnergyError;
use ff_energy::NucleotideVec;

use crate::Move;
use crate::Neighbors;

/// Why a move set cannot be built for a structure.
#[derive(Debug)]
pub enum MoveSetError {
    /// Kinetics support one strand without pseudoknots, linear or circular
    /// (at full length).
    Unsupported,
    /// The energy model cannot evaluate a loop.
    Energy(EnergyError),
}

impl std::fmt::Display for MoveSetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unsupported => write!(f, "Kinetics support a single strand without pseudoknots, a circular strand at full length."),
            Self::Energy(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for MoveSetError {}

/// The moves of a structure, and their changes as moves are applied: what a
/// stochastic simulation needs from a move set.
pub trait MoveSet {
    /// The current structure.
    fn current_structure(&self) -> DotBracketVec;

    fn current_energy(&self) -> i32;

    /// The sequence, with the edits made so far (see `MoveTypes::a_to_i`).
    fn current_sequence(&self) -> &NucleotideVec;

    /// A function to list all possible moves.
    fn propose_moves(&self) -> impl Iterator<Item = (Move, i32)> + '_;

    /// Apply a move, and return the changes to the proposed moves as
    /// (old, new): `old` are moves that no longer exist, `new` are moves
    /// that are new or may have a new energy (unchanged ones included).
    /// A move in `new` that is also proposed already replaces that one.
    fn apply_move(&mut self, mv: &Move) -> (Neighbors, Neighbors);
}

/// A move set on a strand that grows at its 3' end, one base at a time, as
/// during transcription.
pub trait Elongation: MoveSet {
    /// The length of the full sequence.
    fn sequence_length(&self) -> usize;

    /// The length of the current structure.
    fn current_structure_length(&self) -> usize;

    /// Extend the structure by one unpaired base at its 3' end, with the
    /// changes as for `MoveSet::apply_move`.
    fn apply_extension(&mut self) -> (Neighbors, Neighbors);
}
