//! The rates of the moves of the current structure, for drawing the next
//! move of a stochastic simulation in logarithmic time.
//!
//! The moves are the nodes of a complete binary tree, stored as an array
//! (1-based: the children of node i are 2i and 2i + 1). Each node holds its
//! rate and the sum of the rates in its subtree. Drawing a move walks from
//! the root down; changing, adding or removing a move updates the sums on
//! its path to the root. Both take O(log n) steps for n moves. A hash map
//! gives the node of a move.
//!
//! The sums are computed again from the children at every update, never
//! by adding and subtracting differences, so rounding errors do not
//! accumulate over a simulation.

use rustc_hash::FxHashMap;

use crate::Move;

#[derive(Clone, Copy, Debug)]
struct Node {
    rate: f64,
    /// The rate of this node and of all nodes below it.
    rate_sum: f64,
    mv: Move,
}

/// The moves with a positive rate, and the sum of the rates of each subtree.
pub struct RateTree {
    /// The nodes, from index 1; index 0 is unused.
    nodes: Vec<Node>,
    /// The index of each move.
    index: FxHashMap<Move, usize>,
}

impl Default for RateTree {
    fn default() -> Self {
        let unused = Node { rate: 0.0, rate_sum: 0.0, mv: Move::Add { i: 0, j: 0 } };
        Self { nodes: vec![unused], index: FxHashMap::default() }
    }
}

/// A tree of the given moves and rates, in O(n).
impl FromIterator<(Move, f64)> for RateTree {
    fn from_iter<I: IntoIterator<Item = (Move, f64)>>(moves: I) -> Self {
        let mut tree = RateTree::default();
        for (mv, rate) in moves {
            debug_assert!(rate > 0.0);
            tree.index.insert(mv, tree.nodes.len());
            tree.nodes.push(Node { rate, rate_sum: rate, mv });
        }
        for i in (1..tree.nodes.len() / 2 + 1).rev() {
            tree.nodes[i].rate_sum = tree.subtree_sum(i);
        }
        tree
    }
}

impl RateTree {
    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.nodes.len() - 1
    }

    #[cfg(test)]
    pub fn is_empty(&self) -> bool {
        self.nodes.len() == 1
    }

    /// The sum of all rates (the flux), 0 without moves.
    pub fn total_rate(&self) -> f64 {
        self.nodes.get(1).map_or(0.0, |n| n.rate_sum)
    }

    #[cfg(test)]
    /// The moves and their rates, in no particular order.
    pub fn iter(&self) -> impl Iterator<Item = (Move, f64)> + '_ {
        self.nodes[1..].iter().map(|n| (n.mv, n.rate))
    }

    /// Add a move that is not in the tree.
    pub fn insert(&mut self, mv: Move, rate: f64) {
        debug_assert!(rate > 0.0);
        debug_assert!(!self.index.contains_key(&mv));
        let i = self.nodes.len();
        self.nodes.push(Node { rate, rate_sum: rate, mv });
        self.index.insert(mv, i);
        self.update_sums(i / 2);
    }

    /// Change the rate of a move; false if the move is not in the tree.
    pub fn update_rate(&mut self, mv: &Move, rate: f64) -> bool {
        debug_assert!(rate > 0.0);
        let Some(&i) = self.index.get(mv) else {
            return false;
        };
        self.nodes[i].rate = rate;
        self.update_sums(i);
        true
    }

    /// Put a new move with its rate in the place of an old one, without
    /// changing the shape of the tree; false if the old move is not in it.
    pub fn replace(&mut self, old: &Move, new: &Move, rate: f64) -> bool {
        debug_assert!(rate > 0.0);
        let Some(i) = self.index.remove(old) else {
            return false;
        };
        debug_assert!(!self.index.contains_key(new));
        self.index.insert(*new, i);
        self.nodes[i].mv = *new;
        self.nodes[i].rate = rate;
        self.update_sums(i);
        true
    }

    /// Remove a move: the last node takes its place. False if the move is
    /// not in the tree.
    pub fn remove(&mut self, mv: &Move) -> bool {
        let Some(i) = self.index.remove(mv) else {
            return false;
        };
        let last = self.nodes.len() - 1;
        if i != last {
            self.nodes[i] = self.nodes[last];
            self.index.insert(self.nodes[i].mv, i);
            // Includes the last node twice until it is removed below; its
            // path to the root is updated then.
            self.update_sums(i);
        }
        self.nodes.pop();
        self.update_sums(last / 2);
        true
    }

    /// The move at a threshold in [0, total_rate): the rates are laid out
    /// in the order node, left subtree, right subtree.
    pub fn select_by_threshold(&self, mut threshold: f64) -> Option<Move> {
        let mut i = 1;
        while i < self.nodes.len() {
            let node = &self.nodes[i];
            threshold -= node.rate;
            if threshold <= 0.0 {
                return Some(node.mv);
            }
            let left = 2 * i;
            if left >= self.nodes.len() {
                break;
            }
            if threshold < self.nodes[left].rate_sum {
                i = left;
            } else {
                threshold -= self.nodes[left].rate_sum;
                i = left + 1;
            }
        }
        // Only reached if rounding makes the threshold exceed the sums.
        log::warn!("RateTree: threshold beyond the total rate (rounding)");
        self.nodes[1..].last().map(|n| n.mv)
    }

    /// The rate of node i plus the sums of its children.
    fn subtree_sum(&self, i: usize) -> f64 {
        let mut sum = self.nodes[i].rate;
        if let Some(left) = self.nodes.get(2 * i) {
            sum += left.rate_sum;
            if let Some(right) = self.nodes.get(2 * i + 1) {
                sum += right.rate_sum;
            }
        }
        sum
    }

    /// Compute the sums again from node i up to the root. Stops early at a
    /// node with children whose sum did not change: the sums above it are
    /// computed from the same values.
    fn update_sums(&mut self, mut i: usize) {
        while i >= 1 {
            let sum = self.subtree_sum(i);
            if 2 * i < self.nodes.len() && self.nodes[i].rate_sum == sum {
                break;
            }
            self.nodes[i].rate_sum = sum;
            i /= 2;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{RngExt, SeedableRng};
    use rand::rngs::StdRng;

    /// Every sum is that of its subtree, the index points to each node, and
    /// the move at the middle of each interval (node, left, right) is the
    /// one of that interval.
    fn check(tree: &RateTree, model: &FxHashMap<Move, f64>) {
        assert_eq!(tree.len(), model.len());
        for i in 1..tree.nodes.len() {
            assert_eq!(tree.nodes[i].rate_sum, tree.subtree_sum(i));
            assert_eq!(tree.index[&tree.nodes[i].mv], i);
            assert_eq!(model[&tree.nodes[i].mv], tree.nodes[i].rate);
        }
        fn preorder(tree: &RateTree, i: usize, start: &mut f64) {
            if i >= tree.nodes.len() {
                return;
            }
            let rate = tree.nodes[i].rate;
            let mid = *start + rate / 2.0;
            assert_eq!(tree.select_by_threshold(mid), Some(tree.nodes[i].mv));
            *start += rate;
            preorder(tree, 2 * i, start);
            preorder(tree, 2 * i + 1, start);
        }
        let mut start = 0.0;
        preorder(tree, 1, &mut start);
    }

    #[test]
    fn random_operations_as_a_map() {
        let mut rng = StdRng::seed_from_u64(1);
        let mv = |k: usize| Move::Add { i: k, j: k + 4 };
        let initial: Vec<(Move, f64)> = (0..20).map(|k| (mv(k), rng.random_range(0.1..10.0))).collect();
        let mut model: FxHashMap<Move, f64> = initial.iter().copied().collect();
        let mut tree: RateTree = initial.into_iter().collect();
        check(&tree, &model);

        for _ in 0..2000 {
            let k = rng.random_range(0..40);
            let rate = rng.random_range(1e-3..10.0);
            match rng.random_range(0..4) {
                0 => {
                    if model.insert(mv(k), rate).is_none() {
                        tree.insert(mv(k), rate);
                    } else {
                        assert!(tree.update_rate(&mv(k), rate));
                    }
                }
                1 => {
                    let present = model.contains_key(&mv(k));
                    if present {
                        model.insert(mv(k), rate);
                    }
                    assert_eq!(tree.update_rate(&mv(k), rate), present);
                }
                2 => assert_eq!(tree.remove(&mv(k)), model.remove(&mv(k)).is_some()),
                _ => {
                    let new = mv(k + 40);
                    if !model.contains_key(&new) {
                        let present = model.remove(&mv(k)).is_some();
                        if present {
                            model.insert(new, rate);
                        }
                        assert_eq!(tree.replace(&mv(k), &new, rate), present);
                    }
                }
            }
            check(&tree, &model);
            let total: f64 = model.values().sum();
            assert!((tree.total_rate() - total).abs() <= 1e-9 * total.max(1.0));
        }
    }
}
