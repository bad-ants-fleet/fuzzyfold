
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::Move;
use crate::AddDelMoves;

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

impl<'a, E: EnergyModel> Walker for AddDelMoves<'a, E> {
    fn len(&self) -> usize {
        self.del_neighbors().len() + self.add_neighbors().values().map(|v| v.len()).sum::<usize>()
    }

    fn is_empty(&self) -> bool {
        self.del_neighbors().is_empty() && self.add_neighbors().is_empty()
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
            .flat_map(|moves| moves.iter())
            .map(|&(i, j, d)| (Move::Add { i, j }, d)
            );
        let del_moves = self.del_neighbors()
            .iter()
            .map(move |(&i, &delta_e)| {
                (Move::Del { i, j: ltab.pair_lookup(&i) }, delta_e)
            });
        add_moves.chain(del_moves)
    }

    fn apply_move(&mut self, mv: &Move) -> (Moves, Moves) {
        match &mv {
            Move::Add { i, j } => { 
                self.apply_add_move(*i, *j)
            },
            Move::Del { i, j } => { 
                self.apply_del_move(*i, *j)
            },
            _ => panic!("wrong move type!"),
        }
    }
}

