use nohash_hasher::IntMap;

use ff_energy::NestedLoop;
use ff_energy::EnergyModel;

use crate::Move;
use crate::Neighbors;
use crate::movesets::loop_table::LoopTable;

type Pair = (usize, usize);


/// Three-way shift move enumeration.
///
/// Note, this code relies on the existance of a LoopTable, which maps the
/// IntMap indices to the actual NestedLoop.
///
/// Three-way shift moves turn two existing loops into two new loops: 
/// (outer0, inner0) -> (outer1, inner1). The initial pair (i, j) is
/// used to determine which loop is the outer "enclosing" loop, and 
/// which one is the inner "enclosed" loop.
///
/// The activation energy of a three-way shift move is calculated
/// as the maximum over three potential energy barries:
///  -- (E(inner1) + E(outer1)) - (E(inner0) + E(outer0))
///  -- E(inner1) - E(inner0)
///  -- E(outer1) - E(outer0)
///
/// The two loops on either side of a pair, (outer, inner), and the loop
/// without the pair: a three-way shift of the pair splits it again.
struct Joined {
    pair: Pair,
    outer: usize,
    inner: usize,
    combo: NestedLoop,
}

impl Joined {
    fn of<E: EnergyModel>(ltab: &LoopTable<E>, (i, j): Pair) -> Self {
        let (outer, inner) = (ltab.loop_lookup(i), ltab.loop_lookup(j));
        let combo = ltab.get(outer).0.join_loop(&ltab.get(inner).0);
        Joined { pair: (i, j), outer, inner, combo }
    }
}

#[derive(Clone, Default)]
pub struct ThreeWayNeighbors {
    map: IntMap<usize, Neighbors>,
} 

impl ThreeWayNeighbors {
    pub fn map(&self) -> &IntMap<usize, Neighbors> {
        &self.map
    }

    pub fn remove(&mut self, index: &usize) -> Option<Neighbors> {
        self.map.remove(index)
    }

    // Calculate neighbors and return reference.
    pub fn compute_neighbors<E: EnergyModel>(
        &mut self,
        ltab: &LoopTable<E>,
        index: usize,
    ) -> &Neighbors {
        let (init, _) = ltab.get(index);
        let pairs = init.index_pairs();

        let mut loopdict: IntMap<usize, Pair> = IntMap::default();
        // Each pair joined out once, for all of its shifts.
        let mut joined: IntMap<usize, Joined> = IntMap::default();
        for (i, j) in pairs.into_iter() {
            loopdict.insert(i, (i, j));
            loopdict.insert(j, (i, j));
            joined.insert(i, Joined::of(ltab, (i, j)));
        }

        let mut neighbors: Neighbors = Vec::new();
        for s in init.stretches() {
            let ((p5, p3), unpaired) = (s.ends, s.unpaired);
            self.shift_iter(ltab, (p5, unpaired.start), (p3, unpaired.end), &loopdict, &joined, &mut neighbors);
        }
        self.map.insert(index, neighbors);
        &self.map[&index]
    }
    
    #[inline]
    fn shift_iter<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        (p5, p5s): (usize, usize),
        (p3, p3e): (usize, usize),
        loopdict: &IntMap<usize, Pair>,
        joined: &IntMap<usize, Joined>,
        neighbors: &mut Neighbors,
    ) {
        for k in p5s..p3e {
            for (&p, &(i, j)) in loopdict.iter() {
                if !ltab.can_pair(k, p) {
                    continue;
                }
                if (p == p5 && p + ltab.min_hairpin_size() >= k) 
                    || (p == p3 && k + ltab.min_hairpin_size() >= p) 
                {
                    continue;
                } 
                let mv = if p == i {
                    Move::ShiftIK { i, j, k }
                } else {
                    Move::ShiftJK { i, j, k }
                };
                if let Some(delta) = self.get_activation_energy(ltab, &joined[&i], &mv) {
                    neighbors.push((mv, delta))
                }
            }
        }
    }
 
    /// None if a new loop does not fit (see `LoopTable::fits`).
    fn get_activation_energy<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        joined: &Joined,
        mv: &Move,
    ) -> Option<i32> {
        debug_assert_eq!(joined.pair, mv.deleted_pair());
        let (outer0_idx, inner0_idx) = (joined.outer, joined.inner);
        let (p, q) = mv.added_pair();
        let (outer1, inner1) = joined.combo.split_loop(p, q);
        if !ltab.fits(&outer1) || !ltab.fits(&inner1) {
            return None;
        }

        let (_o, outer0_en) = ltab.get(outer0_idx);
        let (_i, inner0_en) = ltab.get(inner0_idx);
        let outer1_en = ltab.energy_of_loop(&outer1);
        let inner1_en = ltab.energy_of_loop(&inner1);

        #[cfg(feature = "shift_analysis")] // For a Figure in Badelt et al. 2026
        {
            let d_en = (outer1_en + inner1_en) - (outer0_en + inner0_en);
            let a_en = (outer1_en - outer0_en).max(inner1_en - inner0_en);
            let bm = match (_o.is_stack(), _i.is_stack(), outer1.is_stack(), inner1.is_stack()) 
            {
                (true, false, false, true) => "neutral_trans",
                (false, true, true, false) => "neutral_trans",
                (false, true, false, true) => "neutral_cis", // loop
                (true, false, false, false) => "del_stack",
                (false, true, false, false) => "del_stack",
                (false, false, true, false) => "add_stack",
                (false, false, false, true) => "add_stack",
                (false, false, false, false) => "other",
                _ => panic!("strange rearrangement"),
            };
            println!("tw {} {} {} {}", d_en, a_en, 0.max(d_en).max(a_en), bm);
        }

        Some(((outer1_en + inner1_en) - (outer0_en + inner0_en))
            .max(outer1_en - outer0_en)
            .max(inner1_en - inner0_en))
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use ff_structure::MultiPairTable;
    use ff_energy::Fuzzyfold;
    use ff_energy::NucleotideVec;

    macro_rules! setup_loop_table {
        ($name:ident, $seq:expr, $db:expr) => {
            let model = Fuzzyfold::default();
            let sequence = NucleotideVec::try_from($seq).unwrap();
            let pairings = MultiPairTable::try_from($db)
                .expect("Invalid structure");
            let $name = LoopTable::try_from((sequence, &pairings, Arc::new(model)))
                .expect("Invalid sequence/structure combination");
        };
    }

    #[test]
    fn test_three_way_setup() {
        setup_loop_table!(ltab, "GUACACGUCG", 
                                "..........");
        let mut twn = ThreeWayNeighbors::default();
        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.is_empty(), "Expected no three-way shift neighbors.");
    }
 
    #[test]
    fn test_three_way_simple() {
        setup_loop_table!(ltab, "GUACCCCUCG", 
                                "...(.....)");
        let mut twn = ThreeWayNeighbors::default();

        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.len() == 2);
        let m1 = Move::ShiftJK { i: 3, j: 9, k: 4 };
        let m2 = Move::ShiftJK { i: 3, j: 9, k: 5 };
        assert_eq!(res[0].0, m1);
        assert_eq!(res[1].0, m2);

        let res = twn.compute_neighbors(&ltab, 1);
        assert!(res.len() == 1);
        let m = Move::ShiftJK { i: 3, j: 9, k: 1 };
        assert_eq!(res[0].0, m);
    }
  
    #[test]
    fn test_three_way_mini() {
        setup_loop_table!(ltab, "CCCCUCG", 
                                "(.....)");
        let mut twn = ThreeWayNeighbors::default();

        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.len() == 2);
        let m1 = Move::ShiftJK { i: 0, j: 6, k: 1 };
        let m2 = Move::ShiftJK { i: 0, j: 6, k: 2 };
        assert_eq!(res[0].0, m1);
        assert_eq!(res[1].0, m2);

        let res = twn.compute_neighbors(&ltab, 1);
        assert!(res.is_empty());
    }
}
