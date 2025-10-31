use std::collections::VecDeque;
use std::usize;
use nohash_hasher::IntMap;
use rustc_hash::{FxHashMap, FxHashSet};

use ff_structure::{LoopInfo, LoopTable, PairTable};
use ff_structure::NAIDX;
use crate::{PairList, Pair};

#[derive(Debug, Clone)]
pub struct PartialOrder {
    pairs: FxHashSet<Pair>,
    precedence: FxHashMap<Pair, FxHashSet<Pair>>,   // DAG: a -> b means a < b (b is a successor)
    predecessors: FxHashMap<Pair, FxHashSet<Pair>>, // DAG: a -> b means b < a (a is a predecessor
    pair_tables: IntMap<usize, PairTable>, // level -> pair_table
}

impl PartialOrder {
    pub fn new() -> Self {
        Self {
            pairs: FxHashSet::default(),
            precedence: FxHashMap::default(),
            predecessors: FxHashMap::default(),
            pair_tables: IntMap::default(),
        }
    }

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
                    return true;
                } else {
                    eprintln!("Warning: missing previous pair table of length {}", length - 1);
                    return false;
                }
            }
        };

        let new_pairs = PairList::try_from(pair_table).expect("PairTable must be valid").pairs;
        for &pair in &new_pairs {
            let _ = self.pairs.insert(pair);
        }
        
        // Now make sure non of the pairs can change anything in the history of the path.
        for (&len, pt) in &self.pair_tables {
            for &pair in &new_pairs {
                if pair.1 as usize > len {
                    continue
                }
                let mut copy = pt.clone();
                match apply_pair_to_pt(&mut copy, pair, &self.precedence) {
                    Ok(Some(old)) => {
                        if old == pair {
                            continue;
                        }
                        // pair < old! Otherwise it would mess up earlier tables!
                        self.precedence.entry(pair).or_default().insert(old); 
                        self.predecessors.entry(old).or_default().insert(pair);
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
        let mut current_pt = extend_pair_table(prev_pt);
        if !self.resolve_conflicts(&new_pairs, &mut current_pt) {
            return false;
        }

        if current_pt.0 != pair_table.0 {
            //eprintln!("Warning: final pair table does not match input");
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
            node: &Pair,
            graph: &FxHashMap<Pair, FxHashSet<Pair>>,
            visited: &mut FxHashSet<Pair>,
            stack: &mut FxHashSet<Pair>,
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

        let mut visited = FxHashSet::default();
        let mut stack = FxHashSet::default();

        for node in self.pairs.iter() {
            if find_cycle_dfs(node, &self.precedence, &mut visited, &mut stack) {
                return false;
            }
        }
        true
    }

    fn resolve_conflicts(&mut self, pairs: &[Pair], pt: &mut PairTable,
    ) -> bool {
        let mut queue: VecDeque<Pair> = pairs.iter().copied().rev().collect();
        let mut progress = true;

        while progress && !queue.is_empty() {
            progress = false;
            let mut skipped = VecDeque::new();

            while let Some(pair) = queue.pop_front() {
                // If it applies, it must be save to apply!
                match apply_pair_to_pt(pt, pair, &self.precedence) {
                    Ok(Some(old)) => {
                        progress = true;
                        if old == pair {
                            continue;
                        }
                        // old < pair! We are now save to apply the move.
                        self.precedence.entry(old).or_default().insert(pair);
                        self.predecessors.entry(pair).or_default().insert(old);
                    }
                    Ok(None) => {
                        progress = true;
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

    pub fn pair_hierarchy(&self) -> FxHashMap<(NAIDX, NAIDX), usize> {
        // Pairs with no predecessors are roots
        let mut levels: FxHashMap<(NAIDX, NAIDX), usize> = FxHashMap::default();
        let mut queue: VecDeque<Pair> = self.pairs.iter()
            .filter(|e| !self.predecessors.contains_key(e))
            .copied()
            .collect();

        for &root in &queue {
            levels.insert((root.0, root.1), 1);
        }

        let mut debug: usize = 0;
        while let Some(edge) = queue.pop_front() {
            let level = levels[&(edge.0, edge.1)];
            if let Some(children) = self.precedence.get(&edge) {
                for &child in children {
                    let key = (child.0, child.1);
                    let child_level = levels.get(&key).copied().unwrap_or(0);

                    if level + 1 > child_level {
                        levels.insert(key, level + 1);
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

    pub fn all_total_orders(&self) -> Vec<Vec<Pair>> {
        let mut all = Vec::new();
        let mut current = Vec::new();
        let mut in_deg: FxHashMap<Pair, usize> = FxHashMap::default();

        for e in &self.pairs {
            in_deg.entry(*e).or_insert(0);
        }

        for (_, targets) in &self.precedence {
            for tgt in targets {
                *in_deg.entry(*tgt).or_insert(0) += 1;
            }
        }

        let mut available: FxHashSet<Pair> = in_deg
            .iter()
            .filter_map(|(&e, &deg)| if deg == 0 { Some(e) } else { None })
            .collect();

        Self::dfs(&self.precedence, &mut in_deg, &mut available, &mut current, &mut all);
        all
    }

    fn dfs(
        graph: &FxHashMap<Pair, FxHashSet<Pair>>,
        in_deg: &mut FxHashMap<Pair, usize>,
        available: &mut FxHashSet<Pair>,
        current: &mut Vec<Pair>,
        all: &mut Vec<Vec<Pair>>,
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

fn extend_pair_table(prev: &PairTable) -> PairTable {
    let mut pt = prev.clone();
    pt.0.push(None);
    pt
}

fn apply_pair_to_pt(pt: &mut PairTable, pair: Pair, pred: &FxHashMap<Pair, FxHashSet<Pair>>,
) -> Result<Option<Pair>, String> {
    use LoopInfo::*;

    let (i, j) = (pair.0 as usize - 1, pair.1 as usize- 1);
    assert!(i<j);

    if Some(j as NAIDX) == pt[i] && Some(i as NAIDX) == pt[j] {
        return Ok(Some(pair));
    }

    let lt = LoopTable::try_from(&*pt).expect("Must be a valid PairTable");
    match (lt[i], lt[j]) {
        (Unpaired { l: iloop }, Unpaired { l: jloop }) => {
            if iloop == jloop {
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                return Ok(None);
            } else {
                Err("Unpaired bases are in different loops.".to_string())
            }
        }
        (Unpaired { l: iloop }, Paired { i: inner_loop, o: outer_loop }) => {
            if iloop == inner_loop || iloop == outer_loop {
                let pi = pt[j].unwrap() as usize;
                let p_edge = if pi < j { Pair(pi as NAIDX+1, j as NAIDX+1) } else { Pair(j as NAIDX+1, pi as NAIDX+1) };
                if pred.get(&pair).map_or(false, |s| s.contains(&p_edge)) {
                    return Err(format!("Precedence violation: ({i}, {j}) < ({pi}, {j})."));
                }
                pt[pi] = None;
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                return Ok(Some(p_edge));
            } else {
                Err(format!("Loop mismatch ({i} unpaired, {j} paired)."))
            }
        }
        (Paired { i: inner_loop, o: outer_loop }, Unpaired { l: jloop }) => {
            if jloop == inner_loop || jloop == outer_loop {
                let pj = pt[i].unwrap() as usize;
                let p_edge = if pj < i { Pair(pj as NAIDX+1, i as NAIDX+1) } else { Pair(i as NAIDX+1, pj as NAIDX+1) };
                if pred.get(&pair).map_or(false, |s| s.contains(&p_edge)) {
                    return Err(format!("Precedence violation: ({i}, {j}) < ({i}, {pj})."));
                }
                pt[pj] = None;
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                return Ok(Some(p_edge));
            } else {
                Err(format!("Loop mismatch ({i} paired, {j} unpaired)."))
            }
        }
        (Paired { i: i_inner, o: i_outer }, Paired { i: j_inner, o: j_outer }) => {
            assert!(i_inner != j_inner);
            if i_outer == j_outer {
                return Err(format!("Both bases paired, but could work. {} {}", i, j));
            } else if i_outer == j_inner {
                return Err(format!("Both bases paired, but could work. {} {}", i, j));
            } else if j_outer == i_inner {
                return Err(format!("Both bases paired, but could work. {} {}", i, j));
            } else {
                return Err(format!("Both bases paired and loop mismatch! {:?}", pair));
            }
        }
        //_  => Err("Cannot form pair: both bases already paired.".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_precedence() {
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        assert!(r);

        println!("{:?}", po.precedence);
        println!("{:?}", po.predecessors);
        assert!(po.precedence.get(&Pair(1,2)).is_none());
        assert!(po.predecessors.get(&Pair(1,2)).is_none());
        assert!(po.precedence.get(&Pair(3,4)).is_none());
        assert!(po.predecessors.get(&Pair(3,4)).is_none());

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1, 2)), Some(&1));
        assert_eq!(ph.get(&(3, 4)), Some(&1));
    }

    #[test]
    fn test_base_precedence_01() {
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
        assert!(r);

        println!("{:?}", po.precedence);
        println!("{:?}", po.predecessors);
        assert!(po.precedence.get(&Pair(1, 2)).unwrap().contains(&Pair(2, 3)));
        assert!(po.predecessors.get(&Pair(2, 3)).unwrap().contains(&Pair(1, 2)));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1,2)), Some(&1));
        assert_eq!(ph.get(&(2,3)), Some(&2));
    }

    #[test]
    fn test_base_precedence_02() {
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        assert!(r);

        println!("{:?}", po.precedence);
        println!("{:?}", po.predecessors);
        assert!(po.precedence.get(&Pair(1, 2)).unwrap().contains(&Pair(1, 3)));
        assert!(po.predecessors.get(&Pair(1, 3)).unwrap().contains(&Pair(1, 2)));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1,2)), Some(&1));
        assert_eq!(ph.get(&(1,3)), Some(&2));
    }

    #[test]
    fn test_invalid_order_01() {
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
        assert!(!r); // no more allowed to apply a move that would have been possible earlier?
        println!("{:?}", po.precedence);
        println!("{:?}", po.predecessors);
        assert!(po.precedence.get(&Pair(1,3)).unwrap().contains(&Pair(1,2)));
        assert!(po.predecessors.get(&Pair(1,2)).unwrap().contains(&Pair(1,3)));
    }

    #[test]
    fn test_invalid_circular_propagation() {
        // ., (), .(), ()(), (()).
        let mut po = PartialOrder::new();
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
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(...)").unwrap());
        assert!(!r); // abusing this test a bit.
        let r = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
        assert!(r);

        let ph = po.pair_hierarchy();
        assert_eq!(ph.get(&(1, 2)), Some(&1));
        assert_eq!(ph.get(&(3, 4)), Some(&1));
        assert_eq!(ph.get(&(1, 5)), Some(&2));

        let e1 = Pair(1,2);
        let e2 = Pair(3,4);
        let e3 = Pair(1,5);

        let orders = po.all_total_orders();
        assert_eq!(orders.len(), 3);
        assert!(orders.contains(&vec![e2, e1, e3]));
        assert!(orders.contains(&vec![e1, e3, e2]));
        assert!(orders.contains(&vec![e1, e2, e3]));
        assert!(!orders.contains(&vec![e2, e3, e1]));
    }

    #[test]
    fn test_precedence_propagation_01() {
        // ., (), ()., ()(), (().)
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("(().)").unwrap());
        assert!(r);

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1, 2)), Some(&3));
        assert_eq!(ph.get(&(3, 4)), Some(&1));
        assert_eq!(ph.get(&(1, 5)), Some(&4));
        assert_eq!(ph.get(&(2, 3)), Some(&2));

        let e1 = Pair(1,2);
        let e2 = Pair(3,4);
        let e3 = Pair(1,5);
        let e4 = Pair(2,3);

        // Confirm transitive dependencies are being tracked
        let p = &po.precedence;
        assert!(p.get(&e1).unwrap().contains(&e3));
        assert!(p.get(&e2).unwrap().contains(&e4));
        assert!(p.get(&e4).unwrap().contains(&e1));
        let q = &po.predecessors;
        assert!(q.get(&e3).unwrap().contains(&e1));
        assert!(q.get(&e4).unwrap().contains(&e2));
        assert!(q.get(&e1).unwrap().contains(&e4));

        let orders = po.all_total_orders();
        assert_eq!(orders.len(), 1);
        assert_eq!(orders[0], vec![e2, e4, e1, e3]);
    }

    #[test]
    fn test_precedence_propagation_02() {
        // . () (.) (.). (.)() ((..))
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)()").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("((..))").unwrap());
        assert!(r);

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1,2)), Some(&1));
        assert_eq!(ph.get(&(1,3)), Some(&2));
        assert_eq!(ph.get(&(1,6)), Some(&3));
        assert_eq!(ph.get(&(4,5)), Some(&1));
        assert_eq!(ph.get(&(2,5)), Some(&2));
    }

    #[test]
    fn test_precedence_propagation_04() {
        // . .. .() ..() (.()) ((()))
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("..").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("..()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("((()))").unwrap());
        assert!(r);

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(2,3)), Some(&1));
        assert_eq!(ph.get(&(3,4)), Some(&2));
        assert_eq!(ph.get(&(1,5)), Some(&2));
        assert_eq!(ph.get(&(1,6)), Some(&3));
        assert_eq!(ph.get(&(2,5)), Some(&1));
    }

    #[test]
    fn test_precedence_propagation_05() {
        // . () (.) ()() ()(). ()(())
        let mut po = PartialOrder::new();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("()(())").unwrap());
        assert!(r);

        println!("{:?}", po.precedence);
        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&(1,2)), Some(&1));
        assert_eq!(ph.get(&(1,3)), Some(&2));
        assert_eq!(ph.get(&(3,4)), Some(&3));
        assert_eq!(ph.get(&(4,5)), Some(&1));
        assert_eq!(ph.get(&(3,6)), Some(&4));
    }

}

