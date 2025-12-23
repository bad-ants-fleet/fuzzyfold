use log::info;
use ahash::AHashSet;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::LoopStructure;

#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
pub enum Move {
    Add {
        i: NAIDX,
        j: NAIDX,
    },
    Del {
        i: NAIDX,
        j: NAIDX,
    },
}

impl Move {
    pub fn inverse(self) -> Self {
        match self {
            Move::Add { i, j } => Move::Del { i, j },
            Move::Del { i, j } => Move::Add { i, j },
        }
    }

    pub fn ij(&self) -> (NAIDX, NAIDX) {
        match self {
            Move::Add { i, j } => (*i, *j),
            Move::Del { i, j } => (*i, *j),
        }
    }
}

pub trait ApplyMove {
    fn apply_move(&mut self, mv: Move);
    fn undo_move(&mut self, mv: Move) {
        self.apply_move(mv.inverse());
    }
}

impl ApplyMove for DotBracketVec {
    fn apply_move(&mut self, mv: Move) {
        match mv {
            Move::Add { i, j } => {
                self[i as usize] = DotBracket::Open;
                self[j as usize] = DotBracket::Close;
            }
            Move::Del { i, j } => {
                self[i as usize] = DotBracket::Unpaired;
                self[j as usize] = DotBracket::Unpaired;
            }
        }
    }
}

impl<'a, E: EnergyModel> ApplyMove for LoopStructure<'a, E> { 
    fn apply_move(&mut self, mv: Move) {
        match mv {
            Move::Add { i, j } => {
                self.apply_add_move(i, j);
            }
            Move::Del { i, j } => {
                self.apply_del_move(i, j);
            }
        }
    }
}

#[derive(Debug)]
struct Frame {
    structure: DotBracketVec,
    depth: usize,
    bp_move: Option<Move>,
    max_delta: i32,
}

impl<'a, E: EnergyModel> LoopStructure<'a, E> {
    pub fn all_moves(&self) -> Vec<(Move, i32)> {
        let mut result = Vec::new();

        for (_, add_neighbors) in self.get_add_neighbors_per_loop().iter() {
            for &(i, j, delta) in add_neighbors {
                result.push((Move::Add { i, j }, delta));
            }
        }
        for (i, j, delta) in self.get_del_neighbors() {
            result.push((Move::Del { i, j }, delta));
        }

        result
    }

    pub fn generate_neighbors<F>(&mut self, 
        maxdelta: i32, 
        maxsteps: usize,
        mut callback: F,
    ) where F: FnMut(&DotBracketVec, i32) {
        let root = DotBracketVec::from(&*self);
        let mut path = Vec::new();

        let mut seen = AHashSet::default();
        seen.insert(root.clone());

        let mut stack = Vec::new();
        stack.push(Frame {
            structure: root,
            depth: 0,
            bp_move: None,
            max_delta: maxdelta,
        });

        while let Some(frame) = stack.pop() {
            let mut db = frame.structure.clone();
            if let Some(bp_move) = frame.bp_move {
                while path.len() >= frame.depth {
                    self.undo_move(path.pop().expect("popping"));
                }
                self.apply_move(bp_move);
                path.push(bp_move);
            }
            debug_assert_eq!(DotBracketVec::from(&*self), db);
            callback(&db, self.energy());

            if path.len() == maxsteps {
                info!("Reached maxsteps during neighbor generation.");
                continue;
            }

            for (bp_move, delta) in self.all_moves() {
                if delta > frame.max_delta {
                    continue;
                }
                db.apply_move(bp_move);
                if seen.insert(db.clone()) {
                    stack.push(Frame {
                        structure: db.clone(),
                        depth: frame.depth + 1,
                        bp_move: Some(bp_move),
                        max_delta: frame.max_delta - delta,
                    });
                }
                db.undo_move(bp_move);
            }
        }
    }
}


