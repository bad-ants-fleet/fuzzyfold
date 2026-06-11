use nohash_hasher::IntMap;

use ff_energy::EnergyModel;
use ff_energy::NearestNeighborLoop;

use crate::Move;
use crate::movesets::loop_table::LoopTable;

type Moves = Vec<(Move, i32)>;

/// Four-way shift move enumeration.
///
/// Note, this code relies on the existance of a LoopTable, which maps the
/// IntMap indices to the actual NearestNeighborLoop.
///
/// Four-way shift moves turn three existing loops into three
/// new loops: (center, merge1, merge2) -> (inside, outer1, outer2) where
/// "inside" will include the two merge loops, and center splits into the outer
/// loops.
///
/// The activation energy for a four-way shift move is calculated
/// as the maximum over three potential energy barries:
///  -- (E(inside) + E(outer1) + E(outer2)) - (E(center) + E(merge1) + E(merg2))
///  -- E(outer1) + E(outer2) - E(center)
///  -- E(inside) - E(merge1) - E(merge2)
///
#[derive(Clone, Default)]
pub struct FourWayNeighbors {
    map: IntMap<usize, Moves>,
} 

impl FourWayNeighbors {
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

    /// Calculate neighbors and return reference.
    pub fn compute_neighbors<E: EnergyModel>(
        &mut self,
        ltab: &LoopTable<E>,
        index: usize,
    ) -> &Moves {
        let (init, _) = ltab.get(index);
        let pairs = init.pairs();

        let mut neighbors = Vec::new(); 
        for (p1, &(i, j)) in pairs.iter().enumerate() {
            let ui = i as usize;
            let uj = j as usize;
            for &(k, l) in pairs[p1 + 1..].iter() {
                let uk = k as usize;
                let ul = l as usize;
                if uj + ltab.min_hairpin_size() < uk 
                    && ltab.can_pair(ui, ul) 
                    && ltab.can_pair(uj, uk) 
                {
                    let mv = Move::ShiftILJK { i, j, k, l };
                    let delta = self.get_activation_energy(ltab, &mv);
                    neighbors.push((mv, delta));
                } else if ltab.can_pair(ui, uk) && ltab.can_pair(ul, uj) 
                    && ui + ltab.min_hairpin_size() < uk 
                    && ul + ltab.min_hairpin_size() < uj
                {
                    let mv = Move::ShiftIKLJ { i, j, k, l };
                    let delta = self.get_activation_energy(ltab, &mv);
                    neighbors.push((mv, delta));
                }
            }
        }
        self.map.insert(index, neighbors);
        &self.map[&index]
    }

    /// Look up / constuct all loops for four-way shift moves.
    pub fn get_loops<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> (usize, // center
          usize, // merge1
          usize, // merge2
          NearestNeighborLoop, // inside
          NearestNeighborLoop, // outer1
          NearestNeighborLoop) // outer2
    {
        match mv {
            Move::ShiftIKLJ { i, j, k, l } => {
                let merge1_idx = ltab.loop_lookup(*i as usize);
                let center_idx = ltab.loop_lookup(*j as usize);
                let merge2_idx = ltab.loop_lookup(*l as usize);
                let (merge1, _energy) = ltab.get(merge1_idx);
                let (center, center_en) = ltab.get(center_idx);
                debug_assert_eq!(center_en, &ltab.geti(*k as usize).1);
                let (merge2, _energy) = ltab.get(merge2_idx);

                let inside = merge1.join_loop(center).join_loop(merge2);
                let (inside, outer1) = inside.split_loop(*i, *k);
                let (inside, outer2) = inside.split_loop(*l, *j);
                (center_idx, merge1_idx, merge2_idx, inside, outer1, outer2)
            },
            Move::ShiftILJK { i, j, k, l } => {

                let center_idx = ltab.loop_lookup(*i as usize);
                let merge1_idx = ltab.loop_lookup(*j as usize);
                let merge2_idx = ltab.loop_lookup(*l as usize);
                let (merge1, _energy) = ltab.get(merge1_idx);
                let (center, center_en) = ltab.get(center_idx);
                debug_assert_eq!(center_en, &ltab.geti(*k as usize).1);
                let (merge2, _energy) = ltab.get(merge2_idx);

                let inside = center.join_loop(merge1).join_loop(merge2);
                let (outer1, inside) = inside.split_loop(*i, *l);
                let (inside, outer2) = inside.split_loop(*j, *k);
                (center_idx, merge1_idx, merge2_idx, inside, outer1, outer2)
            },
            _ => panic!("FourWayNeighbors called with wrong move {:?}", mv),
        }
    }
 
    /// The function to caculate activation energy.
    ///
    /// Note that this function does not necessarily return the free energy
    /// difference!!
    ///
    fn get_activation_energy<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> i32 {
        let (center_idx, merge1_idx, merge2_idx, 
             inside, outer1, outer2) = self.get_loops(ltab, mv);

        let (_, center_en) = ltab.get(center_idx);
        let (_m1, merge1_en) = ltab.get(merge1_idx);
        let (_m2, merge2_en) = ltab.get(merge2_idx);

        let inside_en = ltab.energy_of_loop(&inside);
        let outer1_en = ltab.energy_of_loop(&outer1);
        let outer2_en = ltab.energy_of_loop(&outer2);

        // // A hack for a Figure in Badelt et al. 2026
        // let dE = (inside_en + outer1_en + outer2_en) - (center_en + merge1_en + merge2_en);
        // let aE = (outer1_en + outer2_en - center_en).max(inside_en - merge1_en - merge2_en);
        // let bm = match (_m1.is_stack(), _m2.is_stack(), outer1.is_stack(), outer2.is_stack())
        // {
        //     (true, true, true, true) => "neutral_2s",
        //     (true, false, true, false) => "neutral_1s",
        //     (true, false, false, true) => "neutral_1s",
        //     (false, true, true, false) => "neutral_1s",
        //     (false, true, false, true) => "neutral_1s",

        //     (true, false, true, true) => "add_stack",
        //     (false, true, true, true) => "add_stack",
        //     (false, false, true, false) => "add_stack",
        //     (false, false, false, true) => "add_stack",
        //     (false, false, true, true) => "add_stack",

        //     (true, false, false, false) => "del_stack",
        //     (false, true, false, false) => "del_stack",
        //     (true, true, true, false) => "del_stack",
        //     (true, true, false, true) => "del_stack",
        //     (true, true, false, false) => "del_stack",
        //     _ => "other",
        //     //(false, false, false, false) => "other",
        //     //_ => panic!("strange rearrangement"),
        // };
        // //if bm == "other" {
        // //    println!("{:?} {} {} {} {}", mv, dE, aE, 0.max(dE).max(aE), bm);
        // ////    println!("{:?} {:?} {:?} {:?}", outer1, inner1, _o, _i);
        // //    println!("{} {} {} {}", 
        // //        outer1.is_stack(), outer2.is_stack(), _m1.is_stack(), _m2.is_stack());
        // //}
        // println!("fw {} {} {} {}", dE, aE, 0.max(dE).max(aE), bm);

        ((inside_en + outer1_en + outer2_en) - (center_en + merge1_en + merge2_en))
            .max(outer1_en + outer2_en - center_en)
            .max(inside_en - merge1_en - merge2_en)
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
    fn test_four_way_setup() {
        setup_loop_table!(ltab, "GUACACGUCG", 
                                "..........");
        let mut twn = FourWayNeighbors::default();
        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.is_empty(), "Expected no three-way shift neighbors.");
    }
 
    #[test]
    fn test_four_way_simple() {
        setup_loop_table!(ltab, "AGAAAAACAAAGAAACAA", 
                                ".(.....)...(...)..");
        let mut twn = FourWayNeighbors::default();

        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");

        let res = twn.compute_neighbors(&ltab, 1);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
        
        let res = twn.compute_neighbors(&ltab, 2);
        let m = Move::ShiftILJK { i: 1, j: 7, k: 11, l: 15 };
        assert_eq!(res[0].0, m);

        let (_, _, _, inner, _, _) = twn.get_loops(&ltab, &m);
        assert_eq!(inner, NearestNeighborLoop::Interior { closing: (1, 15), inner: (7, 11) });
    }

    #[test]
    fn test_four_way_simple_rev() {
        setup_loop_table!(ltab, "AGAAAAACAAAGAAACAA", 
                                ".(.....(...)...)..");
        let mut twn = FourWayNeighbors::default();

        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
        
        let res = twn.compute_neighbors(&ltab, 1);
        let m = Move::ShiftIKLJ { i: 1, k: 7, l: 11, j: 15 };
        assert_eq!(res[0].0, m);

        let res = twn.compute_neighbors(&ltab, 2);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
    }
 
}
