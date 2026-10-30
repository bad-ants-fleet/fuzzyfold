//! The neighborhood of a structure: the structures reachable with the moves
//! of a move set within an energy range.

use log::info;
use rustc_hash::FxHashMap;
use rustc_hash::FxHashSet;

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;

use crate::Move;
use crate::MoveSet;

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
                self[i] = DotBracket::OPEN;
                self[j] = DotBracket::CLOSE;
            }
            Move::Edit { .. } => {
                let (i, j) = mv.added_pair();
                self[i] = DotBracket::OPEN;
                self[j] = DotBracket::CLOSE;
            }
            Move::Del { i, j } => {
                self[i] = DotBracket::UNPAIRED;
                self[j] = DotBracket::UNPAIRED;
            }
            Move::ShiftIK { i, j, k } => {
                self[j] = DotBracket::UNPAIRED;
                if i < k {
                    self[i] = DotBracket::OPEN;
                    self[k] = DotBracket::CLOSE;
                } else {
                    self[k] = DotBracket::OPEN;
                    self[i] = DotBracket::CLOSE;
                }
            }
            Move::ShiftJK { i, j, k } => {
                self[i] = DotBracket::UNPAIRED;
                if j < k {
                    self[j] = DotBracket::OPEN;
                    self[k] = DotBracket::CLOSE;
                } else {
                    self[k] = DotBracket::OPEN;
                    self[j] = DotBracket::CLOSE;
                }
            }
            Move::ShiftIKLJ { i, j, k, l } => {
                self[i] = DotBracket::OPEN;
                self[k] = DotBracket::CLOSE;
                self[l] = DotBracket::OPEN;
                self[j] = DotBracket::CLOSE;
            } 
            Move::ShiftILJK { i, j, k, l } => {
                self[i] = DotBracket::OPEN;
                self[j] = DotBracket::OPEN;
                self[k] = DotBracket::CLOSE;
                self[l] = DotBracket::CLOSE;
            }
        }
    }
}


/// A structure to visit: reached with `mv` from the structure of the frame
/// that pushed it, after `depth` moves from the start.
#[derive(Debug)]
struct Frame {
    structure: DotBracketVec,
    depth: usize,
    mv: Option<Move>,
}

/// The structures reachable from the current structure of a move set
/// without exceeding its energy by more than `maxdelta`, also at the
/// barriers of shift moves, in at most `maxsteps` moves. Each structure is
/// passed to the callback once, with its energy. The search is depth
/// first; with a step limit, a structure found again in fewer steps is
/// expanded again. The move set returns to the start.
pub fn enumerate_neighbors<M: MoveSet, F>(moves: &mut M, maxdelta: i32, maxsteps: usize, mut callback: F)
where
    F: FnMut(&DotBracketVec, i32),
{
    let root = moves.current_structure();
    let limit = moves.current_energy().saturating_add(maxdelta);
    // The fewest moves found to each structure, and those passed on.
    let mut depth_of: FxHashMap<DotBracketVec, usize> = FxHashMap::default();
    depth_of.insert(root.clone(), 0);
    let mut reported: FxHashSet<DotBracketVec> = FxHashSet::default();
    let mut path: Vec<Move> = Vec::new();
    let mut stack = vec![Frame { structure: root, depth: 0, mv: None }];

    while let Some(frame) = stack.pop() {
        // Found again in fewer moves after this frame was pushed.
        if depth_of[&frame.structure] < frame.depth {
            continue;
        }
        if let Some(mv) = frame.mv {
            while path.len() >= frame.depth {
                let _ = moves.apply_move(&path.pop().expect("popping").inverse());
            }
            let _ = moves.apply_move(&mv);
            path.push(mv);
        }
        let mut db = frame.structure;
        debug_assert_eq!(moves.current_structure(), db);
        if reported.insert(db.clone()) {
            callback(&db, moves.current_energy());
        }
        if frame.depth == maxsteps {
            info!("Reached maxsteps during neighbor generation.");
            continue;
        }
        let budget = limit - moves.current_energy();
        let depth = frame.depth + 1;
        for (mv, delta) in moves.propose_moves().collect::<Vec<_>>() {
            if delta > budget {
                continue;
            }
            db.apply_move(&mv);
            let shorter = match depth_of.get(&db) {
                None => true,
                Some(&d) => maxsteps != usize::MAX && depth < d,
            };
            if shorter {
                depth_of.insert(db.clone(), depth);
                stack.push(Frame { structure: db.clone(), depth, mv: Some(mv) });
            }
            db.undo_move(&mv);
        }
    }
    while let Some(mv) = path.pop() {
        let _ = moves.apply_move(&mv.inverse());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::collections::VecDeque;
    use ff_structure::MultiPairTable;
    use ff_energy::{Fuzzyfold, NucleotideVec};
    use crate::{LoopNeighbors, MoveTypes};

    fn move_set(model: &Arc<Fuzzyfold>, seq: &str, db: &DotBracketVec, shifts: MoveTypes) -> LoopNeighbors<Fuzzyfold> {
        let pairs = MultiPairTable::try_from(db).unwrap();
        LoopNeighbors::try_from((NucleotideVec::try_from(seq).unwrap(), &pairs, model.clone(), shifts)).unwrap()
    }

    /// Breadth first, with the moves of each structure built from scratch:
    /// every structure is found in the fewest moves.
    fn reference(model: &Arc<Fuzzyfold>, seq: &str, db: &str, shifts: MoveTypes,
        maxdelta: i32, maxsteps: usize,
    ) -> FxHashSet<DotBracketVec> {
        let root = DotBracketVec::try_from(db).unwrap();
        let limit = move_set(model, seq, &root, shifts).current_energy() + maxdelta;
        let mut seen = FxHashSet::from_iter([root.clone()]);
        let mut queue = VecDeque::from([(root, 0)]);
        while let Some((s, depth)) = queue.pop_front() {
            if depth == maxsteps {
                continue;
            }
            let moves = move_set(model, seq, &s, shifts);
            for (mv, delta) in moves.propose_moves() {
                if moves.current_energy() + delta > limit {
                    continue;
                }
                let mut t = s.clone();
                t.apply_move(&mv);
                if seen.insert(t.clone()) {
                    queue.push_back((t, depth + 1));
                }
            }
        }
        seen
    }

    #[test]
    fn as_the_reference() {
        let model = Arc::new(Fuzzyfold::default());
        let seq = "UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC";
        let db = ".((((....)))).((((........))))...............";
        for shifts in [MoveTypes::NONE, MoveTypes::SHIFTS] {
            for (maxdelta, maxsteps) in [(300, usize::MAX), (600, 3), (i32::MAX / 2, 2), (400, 4)] {
                let mut moves = move_set(&model, seq, &DotBracketVec::try_from(db).unwrap(), shifts);
                let mut found = Vec::new();
                enumerate_neighbors(&mut moves, maxdelta, maxsteps, |s, _| found.push(s.clone()));
                let set: FxHashSet<_> = found.iter().cloned().collect();
                assert_eq!(set.len(), found.len(), "{shifts:?} {maxdelta} {maxsteps}: reported twice");
                assert_eq!(set, reference(&model, seq, db, shifts, maxdelta, maxsteps),
                    "{shifts:?} {maxdelta} {maxsteps}");
                assert_eq!(moves.current_structure().to_string(), db);
            }
        }
    }
}
