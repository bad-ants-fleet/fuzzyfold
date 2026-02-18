use log::info;
use rustc_hash::FxHashSet;

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;

use crate::Move;
use crate::Walker;
use crate::LoopNeighbors;
use crate::shift_policy::ShiftPolicy;

pub trait ApplyMove {
    fn apply_move(&mut self, mv: &Move);
    fn undo_move(&mut self, mv: &Move) {
        self.apply_move(&mv.inverse());
    }
}

impl ApplyMove for DotBracketVec {
    fn apply_move(&mut self, mv: &Move) {
        match *mv {
            Move::Add { i, j } => {
                self[i as usize] = DotBracket::Open;
                self[j as usize] = DotBracket::Close;
            }
            Move::Del { i, j } => {
                self[i as usize] = DotBracket::Unpaired;
                self[j as usize] = DotBracket::Unpaired;
            }
            Move::ShiftIK { i, j, k } => {
                self[j as usize] = DotBracket::Unpaired;
                if i < k {
                    self[i as usize] = DotBracket::Open;
                    self[k as usize] = DotBracket::Close;
                } else {
                    self[k as usize] = DotBracket::Open;
                    self[i as usize] = DotBracket::Close;
                }
            }
            Move::ShiftJK { i, j, k } => {
                self[i as usize] = DotBracket::Unpaired;
                if j < k {
                    self[j as usize] = DotBracket::Open;
                    self[k as usize] = DotBracket::Close;
                } else {
                    self[k as usize] = DotBracket::Open;
                    self[j as usize] = DotBracket::Close;
                }
            }
            Move::ShiftIKLJ { i, j, k, l } => {
                self[i as usize] = DotBracket::Open;
                self[k as usize] = DotBracket::Close;
                self[l as usize] = DotBracket::Open;
                self[j as usize] = DotBracket::Close;
            } 
            Move::ShiftILJK { i, j, k, l } => {
                self[i as usize] = DotBracket::Open;
                self[j as usize] = DotBracket::Open;
                self[k as usize] = DotBracket::Close;
                self[l as usize] = DotBracket::Close;
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

impl<E: EnergyModel, P: ShiftPolicy> LoopNeighbors<E, P> {
    pub fn all_moves(&self) -> Vec<(Move, i32)> {
        self.propose_moves().collect()
    }

    pub fn generate_neighbors<F>(&mut self, 
        maxdelta: i32, 
        maxsteps: usize,
        mut callback: F,
    ) where F: FnMut(&DotBracketVec, i32) {
        let root = self.current_structure();
        let mut path: Vec<Move> = Vec::new();

        let mut seen = FxHashSet::default();
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
                    let _ = self.apply_move(
                        &path.pop().expect("popping").inverse());
                }
                let _ = self.apply_move(&bp_move);
                path.push(bp_move);
            }
            debug_assert_eq!(self.current_structure(), db);
            callback(&db, self.current_energy());

            if path.len() == maxsteps {
                info!("Reached maxsteps during neighbor generation.");
                continue;
            }

            for (bp_move, delta) in self.all_moves() {
                if delta > frame.max_delta {
                    continue;
                }
                db.apply_move(&bp_move);
                if seen.insert(db.clone()) {
                    stack.push(Frame {
                        structure: db.clone(),
                        depth: frame.depth + 1,
                        bp_move: Some(bp_move),
                        max_delta: frame.max_delta - delta,
                    });
                }
                db.undo_move(&bp_move);
            }
        }
    }
}


