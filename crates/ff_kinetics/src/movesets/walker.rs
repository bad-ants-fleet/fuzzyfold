
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::shift_policy::ShiftPolicy;
use crate::Move;
use crate::LoopNeighbors;

type Moves = Vec<(Move, i32)>;

pub trait Walker {

    fn len(&self) -> usize;
    
    fn is_empty(&self) -> bool;

    /// The current structure.
    fn current_structure(&self) -> DotBracketVec;

    fn current_energy(&self) -> i32;

    /// A function to list all possible moves.
    fn propose_moves(&self) -> impl Iterator<Item = (Move, i32)> + '_;

    /// A function to apply a particular move 
    /// -> returns updates to the proposed_moves: Old (outdated) New
    fn apply_move(&mut self, mv: &Move) -> (Moves, Moves);
}

impl<'a, E: EnergyModel, P: ShiftPolicy> Walker for LoopNeighbors<'a, E, P> {
    fn len(&self) -> usize {
        let tl = if P::three_way() { self.three_way_shift_neighbors().len() } else { 0 };
        let fl = if P::four_way() { self.four_way_shift_neighbors().len() } else { 0 };
        self.del_neighbors().len() 
            + self.add_neighbors().values().map(|v| v.len()).sum::<usize>()
            + tl + fl
    }

    fn is_empty(&self) -> bool {
        let te = if P::three_way() { self.three_way_shift_neighbors().is_empty() } else { true };
        let fe = if P::four_way() { self.four_way_shift_neighbors().is_empty() } else { true };
 
        self.del_neighbors().is_empty() 
            && self.add_neighbors().is_empty()
            && te && fe
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
}

