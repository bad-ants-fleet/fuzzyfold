use ahash::AHashMap;

use crate::explore::Move;

#[derive(Clone, Debug)]
struct MoveNode {
    rate: f64,
    rate_sum: f64,
    mv: Move,
}

/// A binary tree storing all possible moves.
pub struct RateTree {
    /// A 1-based vector of Moves.
    entries: Vec<MoveNode>,
    /// To access the index given the Move.
    pos_map: AHashMap<Move, usize>,
}

impl RateTree {

    /// Initialization of a RateTree. 
    /// As capacity, think about how many moves you expect.
    /// Roughly the number of base-pairs in the MFE?
    pub fn new(capacity: usize) -> Self {
        let mut entries = Vec::with_capacity(capacity + 1);
        entries.push(MoveNode {
            rate: 0.0,
            rate_sum: 0.0,
            mv: Move::Add { i: 0, j: 0 },
        }); 
        Self {
            entries,
            pos_map: AHashMap::with_capacity(capacity),
        }
    }

    pub fn len(&self) -> usize {
        self.entries.len() - 1
    }

    pub fn is_empty(&self) -> bool {
        self.entries.len() == 1
    }

    pub fn total_rate(&self) -> f64 {
        if self.is_empty() {
            0.0
        } else {
            self.entries[1].rate_sum
        }
    }
    
    fn parent_idx(&self, i: usize) -> usize {
        i/2
    }

    fn left_child(&self, i: usize) -> Option<(usize, &MoveNode)>{
        let pos = 2 * i;
        if pos < self.entries.len() {
            Some((pos, &self.entries[pos]))
        } else {
            None
        }
    }

    fn right_child(&self, i: usize) -> Option<(usize, &MoveNode)>{
        let pos = 2 * i + 1;
        if pos < self.entries.len() {
            Some((pos, &self.entries[pos]))
        } else {
            None
        }
    }

    fn update_partial_sums(&mut self, mut i: usize) {
        while i >= 1 {
            let osum = self.entries[i].rate_sum;
            let mut sum = self.entries[i].rate;
            if let Some((_, entry)) = self.left_child(i) {
                sum += entry.rate_sum;
                if let Some((_, entry)) = self.right_child(i) {
                    sum += entry.rate_sum;
                }
            }
            if (sum - osum).abs() < f64::EPSILON {
                break
            }
            self.entries[i].rate_sum = sum;
            if i == 1 {
                break;
            }
            i = self.parent_idx(i);
        }
    }

    pub fn insert(&mut self, mv: Move, rate: f64) {
        let idx = self.entries.len();

        self.entries.push(MoveNode {
            rate,
            rate_sum: rate,
            mv,
        });

        self.pos_map.insert(mv, idx);
        self.update_partial_sums(self.parent_idx(idx));
    }

    pub fn update_rate(&mut self, mv: &Move, new_rate: f64) -> bool {
        if let Some(&idx) = self.pos_map.get(mv) {
            self.entries[idx].rate = new_rate;
            self.update_partial_sums(idx);
            true
        } else {
            false
        }
    }

    pub fn dirty_replace(&mut self, old_mv: &Move, new_mv: &Move) -> bool {
        let idx = match self.pos_map.remove(old_mv) {
            Some(i) => i,
            None => return false,
        };
        self.pos_map.insert(*new_mv, idx);
        self.entries[idx].mv = *new_mv;
        true
    }

    pub fn remove(&mut self, mv: Move) -> bool {
        let idx = match self.pos_map.remove(&mv) {
            Some(i) => i,
            None => return false,
        };

        let last = self.entries.len() - 1;

        if idx != last {
            let last_node = self.entries[last].clone();
            self.pos_map.insert(last_node.mv, idx);
            self.entries[idx] = last_node;
            self.update_partial_sums(idx);
        }

        self.entries.pop();
        if last >= 1 {
            self.update_partial_sums(self.parent_idx(last));
        }
        true
    }

    pub fn select_by_threshold(&self, mut thresh: f64) -> Option<Move> {
        let mut i = 1;
        while i < self.entries.len() {
            let node = &self.entries[i];
            thresh -= node.rate;
            if thresh <= 0.0 {
                return Some(node.mv);
            }
            if let Some((l, entry)) = self.left_child(i) {
                if thresh < entry.rate_sum {
                    i = l;
                    continue;
                } else {
                    thresh -= entry.rate_sum;
                    i = l + 1;
                }
            } else {
                break;
            }
        }
        //println!("RateTree: roundoff error! This should be extremely rare!");
        self.entries
            .last()
            .map(|n| n.mv)
    }
}

