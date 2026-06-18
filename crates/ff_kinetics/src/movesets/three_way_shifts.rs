use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::EnergyModel;
use ff_energy::NearestNeighborLoop;

use crate::Move;
use crate::movesets::loop_table::LoopTable;

type Pair = (NAIDX, NAIDX);
type Moves = Vec<(Move, i32)>;


/// Three-way shift move enumeration.
///
/// Note, this code relies on the existance of a LoopTable, which maps the
/// IntMap indices to the actual NearestNeighborLoop.
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
#[derive(Clone, Default)]
pub struct ThreeWayNeighbors {
    map: IntMap<usize, Moves>,
} 

impl ThreeWayNeighbors {
    pub fn len(&self) -> usize {
        self.map.values().map(|v| v.len()).sum::<usize>()
    }

    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    pub fn map(&self) -> &IntMap<usize, Moves> {
        &self.map
    }

    pub fn remove(&mut self, index: &usize) -> Option<Moves> {
        self.map.remove(index)
    }

    // Calculate neighbors and return reference.
    pub fn compute_neighbors<E: EnergyModel>(
        &mut self,
        ltab: &LoopTable<E>,
        index: usize,
    ) -> &Moves {
        let (init, _) = ltab.get(index);
        let pairs = init.pairs();

        let mut loopdict: IntMap<NAIDX, Pair> = IntMap::default();
        for (i, j) in pairs.into_iter() {
            loopdict.insert(i, (i, j));
            loopdict.insert(j, (i, j));
        }

        let mut neighbors: Moves = Vec::new();
        let spans = init.inclusive_loop_ranges();
        match init {
            NearestNeighborLoop::Exterior { .. } => {
                debug_assert!(!spans.is_empty(), "Exterior loop has no spans");
                let n = spans.len() - 1;
                for (e, r) in spans.iter().enumerate() {
                    let p5 = *r.start() as NAIDX;
                    let p5s = if e == 0 { p5 } else { p5 + 1 };
                    let p3 = *r.end() as NAIDX;
                    let p3e = if e == n { p3 + 1 } else { p3 };
                    self.shift_iter(ltab, (p5, p5s), (p3, p3e), &loopdict, &mut neighbors);
               }
            },
            NearestNeighborLoop::Hairpin { .. }
            | NearestNeighborLoop::Interior { .. } 
            | NearestNeighborLoop::Multibranch { .. } => {
                for r in spans {
                    let p5 = *r.start() as NAIDX;
                    let p5s = p5 + 1;
                    let p3 = *r.end() as NAIDX;
                    let p3e = p3;
                    self.shift_iter(ltab, (p5, p5s), (p3, p3e), &loopdict, &mut neighbors);
                }
            },
            _ => todo!("Loop type currently not supported for three-way shifts!"),
        }
        self.map.insert(index, neighbors);
        &self.map[&index]
    }
    
    #[inline]
    fn shift_iter<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        (p5, p5s): (NAIDX, NAIDX),
        (p3, p3e): (NAIDX, NAIDX),
        loopdict: &IntMap<NAIDX, Pair>,
        neighbors: &mut Moves,
    ) {
        for k in p5s..p3e {
            let uk = k as usize;
            for (&p, &(i, j)) in loopdict.iter() {
                let up = p as usize;
                if !ltab.can_pair(uk, up) {
                    continue;
                }
                if (p == p5 && up + ltab.min_hairpin_size() >= uk) 
                    || (p == p3 && uk + ltab.min_hairpin_size() >= up) 
                {
                    continue;
                } 
                let mv = if p == i {
                    Move::ShiftIK { i, j, k }
                } else {
                    Move::ShiftJK { i, j, k }
                };
                let delta = self.get_activation_energy(ltab, &mv);
                neighbors.push((mv, delta))
            }
        }
    }
 
    pub fn get_loops<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> (usize, // outer0
          usize, // inner0
          NearestNeighborLoop, // outer1
          NearestNeighborLoop) // inner1
    {
        let (i, j) = mv.deleted_pair();
        let outer0_idx = ltab.loop_lookup(i as usize);
        let inner0_idx = ltab.loop_lookup(j as usize);
        let (outer0, _en) = ltab.get(outer0_idx);
        let (inner0, _en) = ltab.get(inner0_idx); 
        let combo = outer0.join_loop(inner0);
        let (p, q) = mv.added_pair();
        let (outer1, inner1) = combo.split_loop(p, q);
        (outer0_idx, inner0_idx, outer1, inner1)
    }

    fn get_activation_energy<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> i32 {
        let (outer0_idx, inner0_idx, 
             outer1, inner1) = self.get_loops(ltab, mv);

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

        ((outer1_en + inner1_en) - (outer0_en + inner0_en))
            .max(outer1_en - outer0_en)
            .max(inner1_en - inner0_en)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;

    macro_rules! setup_loop_table {
        ($name:ident, $seq:expr, $db:expr) => {
            let model = ViennaRNA::default();
            let sequence = NucleotideVec::try_from($seq).unwrap();
            let pairings = PairTable::try_from($db)
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
