
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::Move;
use crate::AddDelShiftMoves;

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

impl<'a, E: EnergyModel> Walker for AddDelShiftMoves<'a, E> {
    fn len(&self) -> usize {
        let tl = if let Some(tw) = self.three_way_shift_neighbors() {
            tw.len() } else { 0 };
        let fl = if let Some(fw) = self.four_way_shift_neighbors() {
            fw.len() } else { 0 };
        self.del_neighbors().len() 
            + self.add_neighbors().values().map(|v| v.len()).sum::<usize>()
            + tl + fl
    }

    fn is_empty(&self) -> bool {
        let te = if let Some(tw) = self.three_way_shift_neighbors() {
            tw.is_empty() } else { true };
        let fe = if let Some(fw) = self.four_way_shift_neighbors() {
            fw.is_empty() } else { true };
 
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
            .iter()
            .flat_map(|tw| {
                tw.map()
                    .values()
                    .flat_map(|moves| moves.iter().copied())
            });

        let fw_moves = self.four_way_shift_neighbors()
            .iter()
            .flat_map(|fw| {
                fw.map()
                    .values()
                    .flat_map(|moves| moves.iter().copied())
            });

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

