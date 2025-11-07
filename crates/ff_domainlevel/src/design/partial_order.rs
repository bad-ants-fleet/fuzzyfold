//! Determine the partial order of pair insertions.
//!
//! PartialOrder: Parses an ACFP line by line (PairTable by PairTable) to
//! build the hierarchy of pair strengths. For each new PairTable, the pairs
//! are tested. We need to make sure that 
//!  1) the new pairs cannot affect any structures in the history of the path.
//!  2) the new pairs can transform the previous structure to the new one.
//!

use std::collections::VecDeque;
use nohash_hasher::IntSet;
use nohash_hasher::IntMap;

use ff_structure::P1KEY;
use ff_structure::Pair;
use ff_structure::PairSet;
use ff_structure::PairTable;

use crate::design::apply_move::ApplyMove;

#[derive(Debug, Clone, Default)]
pub struct PartialOrder {
    all_pairs: IntSet<P1KEY>,
    pair_tables: IntMap<usize, PairTable>,      // level -> pair_table
    smaller_than: IntMap<P1KEY, IntSet<P1KEY>>, // DAG: a -> b means a < b (b is a successor)
    greater_than: IntMap<P1KEY, IntSet<P1KEY>>, // DAG: a -> b means b < a (a is a predecessor
}

impl PartialOrder {

    pub fn extend_by_pairtable(&mut self, pair_table: &PairTable) -> bool {
        let length = pair_table.len();

        // Return false with a warning if this length is already seen
        if self.pair_tables.contains_key(&length) {
            eprintln!("Warning: pair_table of length {length} already exists");
            return false;
        }

        // Require previous table unless this is the first one
        let prev_pt = match self.pair_tables.get(&(length - 1)) {
            Some(pt) => pt,
            None => {
                if self.pair_tables.is_empty() {
                    self.pair_tables.insert(length, pair_table.clone());
                    for &pkey in PairSet::from(pair_table).iter_keys() {
                        self.all_pairs.insert(pkey);
                    }
                    return true;
                } else {
                    eprintln!("Warning: missing previous pair table of length {}", length - 1);
                    return false;
                }
            }
        };

        let pset = PairSet::from(pair_table);
        for &pkey in pset.iter_keys() {
            self.all_pairs.insert(pkey);
        }
        let new_pairs: Vec<Pair> = pset.iter().collect();
        
        // Make sure none of the pairs can change anything in the history of the path.
        for (&len, pt) in &self.pair_tables {
            for &pair in &new_pairs {
                if pair.j() as usize >= len {
                    continue
                }
                match pt.try_move(pair) {
                    Ok(Some(old)) => {
                        if old == pair {
                            continue;
                        }
                        // pair < old! Otherwise it would mess up earlier tables!
                        self.smaller_than.entry(pair.key()).or_default().insert(old.key()); 
                        self.greater_than.entry(old.key()).or_default().insert(pair.key());
                    }
                    Ok(None) => {
                        // if a pair would just insert like that earlier, then 
                        // it actually should have. so: nope.
                        return false
                    }  
                    Err(_) => {
                        // If the pair does not apply, it is not a problem here!
                    }
                }
            }
        }
 
        // Build initial pt with length n+1
        let mut current_pt = prev_pt.clone();
        current_pt.extend_once();
        if !self.apply_all_pairs(&mut current_pt, &new_pairs) {
            return false;
        }
        if &current_pt != pair_table {
            return false;
        }
        if !self.dependencies_form_dag() {
            return false;
        }

        self.pair_tables.insert(length, current_pt);
        true
    }
    
    fn dependencies_form_dag(&self) -> bool {
        fn find_cycle_dfs(
            node: &P1KEY,
            graph: &IntMap<P1KEY, IntSet<P1KEY>>,
            visited: &mut IntSet<P1KEY>,
            stack: &mut IntSet<P1KEY>,
        ) -> bool {
            if stack.contains(node) {
                return true; // cycle
            }
            if visited.contains(node) {
                return false; // already explored
            }
            visited.insert(*node);
            stack.insert(*node);
            if let Some(children) = graph.get(node) {
                for child in children {
                    if find_cycle_dfs(child, graph, visited, stack) {
                        return true;
                    }
                }
            }
            stack.remove(node);
            false
        }

        let mut visited = IntSet::default();
        let mut stack = IntSet::default();
        for pkey in self.all_pairs.iter() {
            if find_cycle_dfs(pkey, &self.smaller_than, &mut visited, &mut stack) {
                return false;
            }
        }
        true
    }

    fn apply_all_pairs(&mut self, pt: &mut PairTable, pairs: &[Pair], 
    ) -> bool {
        let mut queue: VecDeque<Pair> = pairs.iter().copied().rev().collect();
        let mut progress = true;

        while progress && !queue.is_empty() {
            progress = false;
            let mut skipped = VecDeque::new();

            while let Some(pair) = queue.pop_front() {
                // If it applies, it must be save to apply!
                match pt.try_move(pair) {
                    Ok(Some(old)) => {
                        if old == pair {
                            progress = true;
                            continue;
                        }
                        if self.smaller_than.get(&pair.key()).is_some_and(|s| s.contains(&old.key())) {
                            // The pair cannot form as it is known to be weaker than old!
                            skipped.push_back(pair);
                            continue;
                        }
                        // old < pair! We are now save to apply the move.
                        progress = true;
                        self.smaller_than.entry(old.key()).or_default().insert(pair.key());
                        self.greater_than.entry(pair.key()).or_default().insert(old.key());
                        pt.apply_move(Some(old), pair);
                    }
                    Ok(None) => {
                        progress = true;
                        pt.apply_move(None, pair);
                    }
                    Err(_) => {
                        skipped.push_back(pair);
                    }
                }
            }
            queue = skipped;
        }

        queue.is_empty()
    }

    pub fn pair_hierarchy(&self) -> IntMap<P1KEY, usize> {
        // Pairs with no predecessors are roots
        let mut levels: IntMap<P1KEY, usize> = IntMap::default();
        let mut queue: VecDeque<P1KEY> = self.all_pairs.iter()
            .filter(|e| !self.greater_than.contains_key(e))
            .copied()
            .collect();

        for &root in &queue {
            levels.insert(root, 1);
        }

        let mut debug: usize = 0;
        while let Some(pkey) = queue.pop_front() {
            let level = levels[&pkey];
            if let Some(children) = self.smaller_than.get(&pkey) {
                for &child in children {
                    let child_level = levels.get(&child).copied().unwrap_or(0);

                    if level + 1 > child_level {
                        levels.insert(child, level + 1);
                        queue.push_back(child);
                    }
                }
            }
            if debug > 1000 {
                panic!("Queue too long — is there a cycle in dependencies?");
            }
            debug += 1;
        }

        levels
    }

    pub fn all_total_orders(&self) -> Vec<Vec<P1KEY>> {
        let mut all = Vec::new();
        let mut current = Vec::new();
        let mut in_deg: IntMap<P1KEY, usize> = IntMap::default();

        for &pkey in &self.all_pairs {
            in_deg.entry(pkey).or_insert(0);
        }

        for targets in self.smaller_than.values() {
            for &tgt in targets {
                *in_deg.entry(tgt).or_insert(0) += 1;
            }
        }

        let mut available: IntSet<P1KEY> = in_deg
            .iter()
            .filter_map(|(&e, &deg)| if deg == 0 { Some(e) } else { None })
            .collect();

        Self::dfs(&self.smaller_than, &mut in_deg, &mut available, &mut current, &mut all);
        all
    }

    fn dfs(
        graph: &IntMap<P1KEY, IntSet<P1KEY>>,
        in_deg: &mut IntMap<P1KEY, usize>,
        available: &mut IntSet<P1KEY>,
        current: &mut Vec<P1KEY>,
        all: &mut Vec<Vec<P1KEY>>,
    ) {
        if available.is_empty() {
            if in_deg.values().all(|&v| v == 0) {
                all.push(current.clone());
            }
            return;
        }

        let options: Vec<_> = available.iter().cloned().collect();

        for edge in options {
            available.remove(&edge);
            current.push(edge);

            let mut modified = Vec::new();
            if let Some(children) = graph.get(&edge) {
                for child in children {
                    if let Some(deg) = in_deg.get_mut(child) {
                        *deg -= 1;
                        if *deg == 0 {
                            available.insert(*child);
                            modified.push(*child);
                        }
                    }
                }
            }

            Self::dfs(graph, in_deg, available, current, all);

            for child in &modified {
                available.remove(child);
            }
            if let Some(children) = graph.get(&edge) {
                for child in children {
                    *in_deg.get_mut(child).unwrap() += 1;
                }
            }
            current.pop();
            available.insert(edge);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_precedence() {
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        assert!(r);

        println!("{:?}", po.smaller_than);
        println!("{:?}", po.greater_than);
        assert!(!po.smaller_than.contains_key(&Pair::new(0,1).key()));
        assert!(!po.greater_than.contains_key(&Pair::new(0,1).key()));
        assert!(!po.smaller_than.contains_key(&Pair::new(2,3).key()));
        assert!(!po.greater_than.contains_key(&Pair::new(2,3).key()));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&Pair::new(0, 1).key()), Some(&1));
        assert_eq!(ph.get(&Pair::new(2, 3).key()), Some(&1));
    }

    #[test]
    fn test_base_precedence_01() {
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
        assert!(r);

        println!("{:?}", po.smaller_than);
        println!("{:?}", po.greater_than);
        assert!(po.smaller_than.get(&Pair::new(0, 1).key()).unwrap().contains(&Pair::new(1, 2).key()));
        assert!(po.greater_than.get(&Pair::new(1, 2).key()).unwrap().contains(&Pair::new(0, 1).key()));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&Pair::new(0,1).key()), Some(&1));
        assert_eq!(ph.get(&Pair::new(1,2).key()), Some(&2));
    }

    #[test]
    fn test_base_precedence_02() {
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        assert!(r);
        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(0, 2).key();

        println!("{:?}", po.smaller_than);
        println!("{:?}", po.greater_than);
        assert!(po.smaller_than.get(&p1).unwrap().contains(&p2));
        assert!(po.greater_than.get(&p2).unwrap().contains(&p1));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&p1), Some(&1));
        assert_eq!(ph.get(&p2), Some(&2));
    }

    #[test]
    fn test_invalid_order_01() {
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
        assert!(!r); // no more allowed to apply a move that would have been possible earlier?
        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(0, 2).key();
        println!("{:?}", po.smaller_than);
        println!("{:?}", po.greater_than);
        assert!(po.smaller_than.get(&p2).unwrap().contains(&p1));
        assert!(po.greater_than.get(&p1).unwrap().contains(&p2));
    }

    #[test]
    fn test_invalid_circular_propagation() {
        // ., (), .(), ()(), (()).
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(()).").unwrap());
        assert!(!r); // would require 4-way migration.
    }

    #[test]
    fn test_multiple_orders() {
        // ., (), ()., ()(), (...)
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(...)").unwrap());
        assert!(!r); // abusing this test a bit.
        let r = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
        assert!(r);

        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(2, 3).key();
        let p3 = Pair::new(0, 4).key();

        let ph = po.pair_hierarchy();
        assert_eq!(ph.get(&p1), Some(&1));
        assert_eq!(ph.get(&p2), Some(&1));
        assert_eq!(ph.get(&p3), Some(&2));

        let orders = po.all_total_orders();
        assert_eq!(orders.len(), 3);
        assert!(orders.contains(&vec![p2, p1, p3]));
        assert!(orders.contains(&vec![p1, p3, p2]));
        assert!(orders.contains(&vec![p1, p2, p3]));
        assert!(!orders.contains(&vec![p2, p3, p1]));
    }

    #[test]
    fn test_precedence_propagation_01() {
        // ., (), ()., ()(), (().)
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(().)").unwrap());
        assert!(r);

        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(2, 3).key();
        let p3 = Pair::new(0, 4).key();
        let p4 = Pair::new(1, 2).key();

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&p1), Some(&3));
        assert_eq!(ph.get(&p2), Some(&1));
        assert_eq!(ph.get(&p3), Some(&4));
        assert_eq!(ph.get(&p4), Some(&2));

        // Confirm transitive dependencies are being tracked
        let p = &po.smaller_than;
        assert!(p.get(&p1).unwrap().contains(&p3));
        assert!(p.get(&p2).unwrap().contains(&p4));
        assert!(p.get(&p4).unwrap().contains(&p1));
        let q = &po.greater_than;
        assert!(q.get(&p3).unwrap().contains(&p1));
        assert!(q.get(&p4).unwrap().contains(&p2));
        assert!(q.get(&p1).unwrap().contains(&p4));

        let orders = po.all_total_orders();
        assert_eq!(orders.len(), 1);
        assert_eq!(orders[0], [p2, p4, p1, p3]);
    }

    #[test]
    fn test_precedence_propagation_02() {
        // . () (.) (.). (.)() ((..))
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("((..))").unwrap());
        assert!(r);

        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(0, 2).key();
        let p3 = Pair::new(0, 5).key();
        let p4 = Pair::new(3, 4).key();
        let p5 = Pair::new(1, 4).key();

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&p1), Some(&1));
        assert_eq!(ph.get(&p2), Some(&2));
        assert_eq!(ph.get(&p3), Some(&3));
        assert_eq!(ph.get(&p4), Some(&1));
        assert_eq!(ph.get(&p5), Some(&2));
    }

    #[test]
    fn test_precedence_propagation_04() {
        // . .. .() ..() (.()) ((()))
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("..").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("..()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("((()))").unwrap());
        assert!(r);

        let p1 = Pair::new(1, 2).key();
        let p2 = Pair::new(2, 3).key();
        let p3 = Pair::new(0, 4).key();
        let p4 = Pair::new(0, 5).key();
        let p5 = Pair::new(1, 4).key();

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&p1), Some(&1));
        assert_eq!(ph.get(&p2), Some(&2));
        assert_eq!(ph.get(&p3), Some(&2));
        assert_eq!(ph.get(&p4), Some(&3));
        assert_eq!(ph.get(&p5), Some(&1));
    }

    #[test]
    fn test_precedence_propagation_05() {
        // . () (.) ()() ()(). ()(())
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("()(())").unwrap());
        assert!(r);

        let p1 = Pair::new(0, 1).key();
        let p2 = Pair::new(0, 2).key();
        let p3 = Pair::new(2, 3).key();
        let p4 = Pair::new(3, 4).key();
        let p5 = Pair::new(2, 5).key();

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&p1), Some(&1));
        assert_eq!(ph.get(&p2), Some(&2));
        assert_eq!(ph.get(&p3), Some(&3));
        assert_eq!(ph.get(&p4), Some(&1));
        assert_eq!(ph.get(&p5), Some(&4));
    }

}

