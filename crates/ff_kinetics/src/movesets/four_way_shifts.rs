use nohash_hasher::IntMap;

use ff_energy::EnergyModel;
use ff_energy::NearestNeighborLoop;

use crate::Move;
use crate::movesets::LoopTable;

type Moves = Vec<(Move, i32)>;

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

    // Calculate neighbors and return reference.
    pub fn compute_neighbors<E: EnergyModel>(
        &mut self,
        ltab: &LoopTable<E>,
        index: usize,
    ) -> &Moves {
        let (init, _) = ltab.get(index);
        let mut neighbors = Vec::new(); 

        let pairs = init.pairs();
        for (p1, (i, j)) in pairs.iter().enumerate() {
            let ui = *i as usize;
            let uj = *j as usize;
            for (k, l) in pairs[p1 + 1..].iter() {
                let uk = *k as usize;
                let ul = *l as usize;
                //TODO switch?
                if ui + ltab.min_hairpin_size() < uk && ltab.can_pair(ui, uk) 
                    && ul + ltab.min_hairpin_size() < uj && ltab.can_pair(ul, uj) {
                        let mv = Move::ShiftIKLJ { i: *i, j: *j, k: *k, l: *l };
                        let delta = self.get_activation_energy(ltab, &mv);
                        neighbors.push((mv, delta));
                } else if uj + ltab.min_hairpin_size() < uk && ltab.can_pair(uj, uk) 
                    && ltab.can_pair(ui, ul) {
                        let mv = Move::ShiftILJK { i: *i, j: *j, k: *k, l: *l };
                        let delta = self.get_activation_energy(ltab, &mv);
                        neighbors.push((mv, delta));
                }
            }
        }
        self.map.insert(index, neighbors);
        &self.map[&index]
    }

    pub fn get_loops<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> (usize, // init
          usize, // merge1
          usize, // merge2
          NearestNeighborLoop, // inner
          NearestNeighborLoop, // outer1
          NearestNeighborLoop) // outer2
    {
        match mv {
            Move::ShiftIKLJ { i, j, k, l } => {
                let m1_idx = ltab.loop_lookup(*i as usize);
                let it_idx = ltab.loop_lookup(*j as usize);
                let m2_idx = ltab.loop_lookup(*l as usize);
                let (m1, _en) = ltab.get(m1_idx);
                let (it, it_en) = ltab.get(it_idx);
                debug_assert_eq!(it_en, &ltab.geti(*k as usize).1);
                let (m2, _en) = ltab.get(m2_idx);

                let inner = m1.join_loop(it).join_loop(m2);
                let (inner, o1) = inner.split_loop(*i, *k);
                let (inner, o2) = inner.split_loop(*l, *j);
                (it_idx, m1_idx, m2_idx, inner, o1, o2)
            },
            Move::ShiftILJK { i, j, k, l } => {
                let it_idx = ltab.loop_lookup(*i as usize);
                let m1_idx = ltab.loop_lookup(*j as usize);
                let m2_idx = ltab.loop_lookup(*l as usize);
                let (it, it_en) = ltab.get(it_idx);
                let (m1, _en) = ltab.get(m1_idx);
                debug_assert_eq!(it_en, &ltab.geti(*k as usize).1);
                let (m2, _en) = ltab.get(m2_idx);
                let inner = it.join_loop(m1).join_loop(m2);
                let (o1, inner) = inner.split_loop(*i, *l);
                let (inner, o2) = inner.split_loop(*j, *k);
                (it_idx, m1_idx, m2_idx, inner, o1, o2)
 
            },
            _ => panic!("wrong move"),
        }

    }
 
    fn get_activation_energy<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> i32 {
        let (it_idx, m1_idx, m2_idx, inner, o1, o2) = self.get_loops(ltab, mv);

        let (_it, it_en) = ltab.get(it_idx);
        let (_m1, m1_en) = ltab.get(m1_idx);
        let (_m2, m2_en) = ltab.get(m2_idx);

        let ir_en = ltab.energy_of_loop(&inner);
        let o1_en = ltab.energy_of_loop(&o1);
        let o2_en = ltab.energy_of_loop(&o2);

        (o1_en + o2_en - it_en)
            .max(ir_en - m1_en - m2_en)
            .max((ir_en + o1_en + o2_en) - (it_en + m1_en + m2_en))
    }

    pub fn remove(
        &mut self, 
        index: &usize, 
    ) -> Option<Moves> {
        self.map.remove(index)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;

    #[test]
    fn test_four_way_setup() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GUACACGUCG");
        let pairings =       PairTable::try_from("..........").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut twn = FourWayNeighbors::default();
        let res = twn.compute_neighbors(&ltab, 0);
        assert!(res.is_empty(), "Expected no three-way shift neighbors.");
    }
 
    #[test]
    fn test_four_way_simple() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("AGAAAAACAAAGAAACAA");
        let pairings =       PairTable::try_from(".(.....)...(...)..").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        println!("LT: {:?}", ltab);
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
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("AGAAAAACAAAGAAACAA");
        let pairings =       PairTable::try_from(".(.....(...)...)..").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        println!("LT: {:?}", ltab);
        let mut twn = FourWayNeighbors::default();

        let res = twn.compute_neighbors(&ltab, 0);
        println!("{:?}", res);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
        
        let res = twn.compute_neighbors(&ltab, 1);
        println!("{:?}", res);
        let m = Move::ShiftIKLJ { i: 1, k: 7, l: 11, j: 15 };
        assert_eq!(res[0].0, m);

        let res = twn.compute_neighbors(&ltab, 2);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
    }
 
}
