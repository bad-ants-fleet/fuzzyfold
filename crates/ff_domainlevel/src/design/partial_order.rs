use std::collections::VecDeque;
use nohash_hasher::IntSet;
use nohash_hasher::IntMap;

use ff_structure::LoopInfo;
use ff_structure::LoopTable;
use ff_structure::PairTable;
use ff_structure::NAIDX;
use ff_structure::P1KEY;
use ff_structure::Pair;
use ff_structure::PairSet;

#[derive(Debug, Clone, Default)]
pub struct PartialOrder {
    all_pairs: IntSet<P1KEY>,
    pair_tables: IntMap<usize, PairTable>,      // level -> pair_table
    greater_than: IntMap<P1KEY, IntSet<P1KEY>>, // DAG: a -> b means a < b (b is a successor)
    smaller_than: IntMap<P1KEY, IntSet<P1KEY>>, // DAG: a -> b means b < a (a is a predecessor
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
                    //NOTE:: no assignment??
                    self.pair_tables.insert(length, pair_table.clone());
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
        
        // Make sure none of the pairs can change anything in the history of the path.
        for (&len, pt) in &self.pair_tables {
            for pair in pset.iter() {
                if pair.j() as usize >= len {
                    continue
                }
                let mut copy = pt.clone();
                match apply_pair_to_pt(&mut copy, pair, &self.greater_than) {
                    Ok(Some(old)) => {
                        if old == pair {
                            continue;
                        }
                        // pair < old! Otherwise it would mess up earlier tables!
                        self.greater_than.entry(pair.key()).or_default().insert(old.key()); 
                        self.smaller_than.entry(old.key()).or_default().insert(pair.key());
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
        current_pt.append_base();
        let new_pairs: Vec<Pair> = pset.iter().collect();
        if !self.resolve_conflicts(&new_pairs, &mut current_pt) {
            return false;
        }

        if &current_pt != pair_table {
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
            if find_cycle_dfs(pkey, &self.greater_than, &mut visited, &mut stack) {
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
                match apply_pair_to_pt(pt, pair, &self.greater_than) {
                    Ok(Some(old)) => {
                        progress = true;
                        if old == pair {
                            continue;
                        }
                        // old < pair! We are now save to apply the move.
                        self.greater_than.entry(old.key()).or_default().insert(pair.key());
                        self.smaller_than.entry(pair.key()).or_default().insert(old.key());
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

    pub fn pair_hierarchy(&self) -> IntMap<P1KEY, usize> {
        // Pairs with no predecessors are roots
        let mut levels: IntMap<P1KEY, usize> = IntMap::default();
        let mut queue: VecDeque<P1KEY> = self.all_pairs.iter()
            .filter(|e| !self.smaller_than.contains_key(e))
            .copied()
            .collect();

        for &root in &queue {
            levels.insert(root, 1);
        }

        let mut debug: usize = 0;
        while let Some(pkey) = queue.pop_front() {
            let level = levels[&pkey];
            if let Some(children) = self.greater_than.get(&pkey) {
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

        for targets in self.greater_than.values() {
            for &tgt in targets {
                *in_deg.entry(tgt).or_insert(0) += 1;
            }
        }

        let mut available: IntSet<P1KEY> = in_deg
            .iter()
            .filter_map(|(&e, &deg)| if deg == 0 { Some(e) } else { None })
            .collect();

        Self::dfs(&self.greater_than, &mut in_deg, &mut available, &mut current, &mut all);
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

fn extend_pair_table(prev: &PairTable) -> PairTable {
    let mut pt = prev.clone();
    pt.append_base();
    pt
}

fn apply_pair_to_pt(
    pt: &mut PairTable, 
    pair: Pair, 
    pred: &IntMap<P1KEY, IntSet<P1KEY>>,
) -> Result<Option<Pair>, String> {
    use LoopInfo::*;

    let (i, j) = (pair.i() as usize, pair.j() as usize);
    assert!(i<j);

    if Some(j as NAIDX) == pt[i] && Some(i as NAIDX) == pt[j] {
        return Ok(Some(pair));
    }

    let lt = LoopTable::from(&*pt);
    match (lt[i], lt[j]) {
        (Unpaired { l: iloop }, Unpaired { l: jloop }) => {
            if iloop == jloop {
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                Ok(None)
            } else {
                Err("Unpaired bases are in different loops.".to_string())
            }
        }
        (Unpaired { l: iloop }, Paired { i: inner_loop, o: outer_loop }) => {
            if iloop == inner_loop || iloop == outer_loop {
                let pi = pt[j].unwrap() as usize;
                let p_edge = if pi < j { Pair::new(pi as NAIDX, j as NAIDX) 
                                } else { Pair::new(j as NAIDX, pi as NAIDX) };
                if pred.get(&pair.key()).map_or(false, |s| s.contains(&p_edge.key())) {
                    return Err(format!("Precedence violation: ({i}, {j}) < ({pi}, {j})."));
                }
                pt[pi] = None;
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                Ok(Some(p_edge))
            } else {
                Err(format!("Loop mismatch ({i} unpaired, {j} paired)."))
            }
        }
        (Paired { i: inner_loop, o: outer_loop }, Unpaired { l: jloop }) => {
            if jloop == inner_loop || jloop == outer_loop {
                let pj = pt[i].unwrap() as usize;
                let p_edge = if pj < i { Pair::new(pj as NAIDX, i as NAIDX) } else { Pair::new(i as NAIDX, pj as NAIDX) };
                if pred.get(&pair.key()).map_or(false, |s| s.contains(&p_edge.key())) {
                    return Err(format!("Precedence violation: ({i}, {j}) < ({i}, {pj})."));
                }
                pt[pj] = None;
                pt[i] = Some(j as NAIDX);
                pt[j] = Some(i as NAIDX);
                Ok(Some(p_edge))
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
        let mut po = PartialOrder::default();
        let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
        let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
        let r = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
        assert!(r);

        println!("{:?}", po.greater_than);
        println!("{:?}", po.smaller_than);
        assert!(!po.greater_than.contains_key(&Pair::new(0,1).key()));
        assert!(!po.smaller_than.contains_key(&Pair::new(0,1).key()));
        assert!(!po.greater_than.contains_key(&Pair::new(2,3).key()));
        assert!(!po.smaller_than.contains_key(&Pair::new(2,3).key()));

        let ph = po.pair_hierarchy();
        println!("{:?}", ph);
        assert_eq!(ph.get(&Pair::new(0, 1).key()), Some(&1));
        assert_eq!(ph.get(&Pair::new(2, 3).key()), Some(&1));
    }

    //#[test]
    //fn test_base_precedence_01() {
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
    //    assert!(r);

    //    println!("{:?}", po.precedence);
    //    println!("{:?}", po.predecessors);
    //    assert!(po.precedence.get(&Pair::new(0, 1)).unwrap().contains(&Pair::new(1, 2)));
    //    assert!(po.predecessors.get(&Pair::new(1, 2)).unwrap().contains(&Pair::new(0, 1)));

    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(0,1)), Some(&1));
    //    assert_eq!(ph.get(&(1,2)), Some(&2));
    //}

    //#[test]
    //fn test_base_precedence_02() {
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
    //    assert!(r);
    //    let p1 = Pair::new(0, 1);
    //    let p2 = Pair::new(0, 2);

    //    println!("{:?}", po.precedence);
    //    println!("{:?}", po.predecessors);
    //    assert!(po.precedence.get(&p1).unwrap().contains(&p2));
    //    assert!(po.predecessors.get(&p2).unwrap().contains(&p1));

    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(0,1)), Some(&1));
    //    assert_eq!(ph.get(&(0,2)), Some(&2));
    //}

    //#[test]
    //fn test_invalid_order_01() {
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
    //    assert!(!r); // no more allowed to apply a move that would have been possible earlier?
    //    let p1 = Pair::new(0, 1);
    //    let p2 = Pair::new(0, 2);
    //    println!("{:?}", po.precedence);
    //    println!("{:?}", po.predecessors);
    //    assert!(po.precedence.get(&p2).unwrap().contains(&p1));
    //    assert!(po.predecessors.get(&p1).unwrap().contains(&p2));
    //}

    //#[test]
    //fn test_invalid_circular_propagation() {
    //    // ., (), .(), ()(), (()).
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(()).").unwrap());
    //    assert!(!r); // would require 4-way migration.
    //}

    //#[test]
    //fn test_multiple_orders() {
    //    // ., (), ()., ()(), (...)
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(...)").unwrap());
    //    assert!(!r); // abusing this test a bit.
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
    //    assert!(r);

    //    let ph = po.pair_hierarchy();
    //    assert_eq!(ph.get(&(0, 1)), Some(&1));
    //    assert_eq!(ph.get(&(2, 3)), Some(&1));
    //    assert_eq!(ph.get(&(0, 4)), Some(&2));

    //    let e1 = Pair::new(0,1);
    //    let e2 = Pair::new(2,3);
    //    let e3 = Pair::new(0,4);

    //    let orders = po.all_total_orders();
    //    assert_eq!(orders.len(), 3);
    //    assert!(orders.contains(&vec![e2, e1, e3]));
    //    assert!(orders.contains(&vec![e1, e3, e2]));
    //    assert!(orders.contains(&vec![e1, e2, e3]));
    //    assert!(!orders.contains(&vec![e2, e3, e1]));
    //}

    //#[test]
    //fn test_precedence_propagation_01() {
    //    // ., (), ()., ()(), (().)
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("().").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("(().)").unwrap());
    //    assert!(r);

    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(0, 1)), Some(&3));
    //    assert_eq!(ph.get(&(2, 3)), Some(&1));
    //    assert_eq!(ph.get(&(0, 4)), Some(&4));
    //    assert_eq!(ph.get(&(1, 2)), Some(&2));

    //    let e1 = Pair::new(0,1);
    //    let e2 = Pair::new(2,3);
    //    let e3 = Pair::new(0,4);
    //    let e4 = Pair::new(1,2);

    //    // Confirm transitive dependencies are being tracked
    //    let p = &po.precedence;
    //    assert!(p.get(&e1).unwrap().contains(&e3));
    //    assert!(p.get(&e2).unwrap().contains(&e4));
    //    assert!(p.get(&e4).unwrap().contains(&e1));
    //    let q = &po.predecessors;
    //    assert!(q.get(&e3).unwrap().contains(&e1));
    //    assert!(q.get(&e4).unwrap().contains(&e2));
    //    assert!(q.get(&e1).unwrap().contains(&e4));

    //    let orders = po.all_total_orders();
    //    assert_eq!(orders.len(), 1);
    //    assert_eq!(orders[0], vec![e2, e4, e1, e3]);
    //}

    //#[test]
    //fn test_precedence_propagation_02() {
    //    // . () (.) (.). (.)() ((..))
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("(.).").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("(.)()").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("((..))").unwrap());
    //    assert!(r);

    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(0,1)), Some(&1));
    //    assert_eq!(ph.get(&(0,2)), Some(&2));
    //    assert_eq!(ph.get(&(0,5)), Some(&3));
    //    assert_eq!(ph.get(&(3,4)), Some(&1));
    //    assert_eq!(ph.get(&(1,4)), Some(&2));
    //}

    //#[test]
    //fn test_precedence_propagation_04() {
    //    // . .. .() ..() (.()) ((()))
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("..").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("..()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("(.())").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("((()))").unwrap());
    //    assert!(r);

    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(1,2)), Some(&1));
    //    assert_eq!(ph.get(&(2,3)), Some(&2));
    //    assert_eq!(ph.get(&(0,4)), Some(&2));
    //    assert_eq!(ph.get(&(0,5)), Some(&3));
    //    assert_eq!(ph.get(&(1,4)), Some(&1));
    //}

    //#[test]
    //fn test_precedence_propagation_05() {
    //    // . () (.) ()() ()(). ()(())
    //    let mut po = PartialOrder::new();
    //    let _ = po.extend_by_pairtable(&PairTable::try_from(".").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("(.)").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()()").unwrap());
    //    let _ = po.extend_by_pairtable(&PairTable::try_from("()().").unwrap());
    //    let r = po.extend_by_pairtable(&PairTable::try_from("()(())").unwrap());
    //    assert!(r);

    //    println!("{:?}", po.precedence);
    //    let ph = po.pair_hierarchy();
    //    println!("{:?}", ph);
    //    assert_eq!(ph.get(&(0,1)), Some(&1));
    //    assert_eq!(ph.get(&(0,2)), Some(&2));
    //    assert_eq!(ph.get(&(2,3)), Some(&3));
    //    assert_eq!(ph.get(&(3,4)), Some(&1));
    //    assert_eq!(ph.get(&(2,5)), Some(&4));
    //}

}

