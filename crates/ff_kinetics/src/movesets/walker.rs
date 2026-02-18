
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::shift_policy::ShiftPolicy;
use crate::Move;
use crate::LoopNeighbors;

pub type Moves = Vec<(Move, i32)>;

pub trait Walker {
    /// The sequence length.
    fn sequence_length(&self) -> usize;

    /// The current structure length.
    fn current_structure_length(&self) -> usize;

    /// The current structure.
    fn current_structure(&self) -> DotBracketVec;

    fn current_energy(&self) -> i32;

    /// A function to list all possible moves.
    fn propose_moves(&self) -> impl Iterator<Item = (Move, i32)> + '_;

    /// A function to apply a particular move 
    /// -> returns updates to the proposed_moves: Old (outdated) New
    fn apply_move(&mut self, mv: &Move) -> (Moves, Moves);

    fn apply_extension(&mut self) -> (Moves, Moves);
}

impl<E: EnergyModel, P: ShiftPolicy> Walker for LoopNeighbors<E, P> {

    fn sequence_length(&self) -> usize {
        self.loop_table().sequence_length()
    }

    fn current_structure_length(&self) -> usize {
        self.loop_table().lookup_len()
    }

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

    fn apply_move(&mut self, mv: &Move) -> (Moves, Moves) {
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

    fn apply_extension(&mut self) -> (Moves, Moves) {
        self.apply_ext_move()
    }
}

