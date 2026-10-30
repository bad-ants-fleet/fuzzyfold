
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_energy::EnergyError;

use crate::Move;
use crate::LoopNeighbors;

/// The neighbors of a structure in a move set: each move with its energy
/// change, or its barrier for shift moves.
pub type Neighbors = Vec<(Move, i32)>;

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

impl<E: EnergyModel> MoveSet for LoopNeighbors<E> {

    fn current_structure(&self) -> DotBracketVec {
        DotBracketVec::from(self.loop_table())
    }

    fn current_energy(&self) -> i32 {
        self.loop_table().energy()
    }

    fn propose_moves(&self) -> impl Iterator<Item = (Move, i32)> + '_ {
        let ltab = self.loop_table();

        let add_moves = self.add_neighbors()
            .values()
            .flat_map(|moves| moves.iter().cloned());

        let del_moves = self.del_neighbors()
            .iter()
            .map(move |(&i, &delta_e)| {
                (Move::Del { i, j: ltab.pair_lookup(&i) }, delta_e)
            });

        let tw_moves = self.three_way_shift_neighbors()
            .map()
            .values()
            .flat_map(|moves| moves.iter().cloned());

        let fw_moves = self.four_way_shift_neighbors() 
            .map()
            .values()
            .flat_map(|moves| moves.iter().cloned());

        add_moves
            .chain(del_moves)
            .chain(tw_moves)
            .chain(fw_moves)
    }

    fn apply_move(&mut self, mv: &Move) -> (Neighbors, Neighbors) {
        match &mv {
            Move::Add { i, j } => { 
                self.apply_add_move(*i, *j)
            },
            Move::Del { i, j } => { 
                self.apply_del_move(*i, *j)
            },
            mv @ Move::ShiftIK { k, .. } => { 
                self.apply_three_way_shift_move(mv, *k)
            },
            mv @ Move::ShiftJK { k, .. } => { 
                self.apply_three_way_shift_move(mv, *k)
            },
            mv @ Move::ShiftIKLJ { .. } => { 
                self.apply_four_way_shift_move(mv)
            },
            mv @ Move::ShiftILJK { .. } => { 
                self.apply_four_way_shift_move(mv)
            },
        }
    }

}

impl<E: EnergyModel> Elongation for LoopNeighbors<E> {
    fn sequence_length(&self) -> usize {
        self.loop_table().sequence_length()
    }

    fn current_structure_length(&self) -> usize {
        self.loop_table().lookup_len()
    }

    fn apply_extension(&mut self) -> (Neighbors, Neighbors) {
        self.apply_ext_move()
    }
}

