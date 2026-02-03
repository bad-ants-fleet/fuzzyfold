use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::EnergyModel;
use ff_energy::NearestNeighborLoop;

use crate::Move;
use crate::movesets::LoopTable;

type Pair = (NAIDX, NAIDX);
type Moves = Vec<(Move, i32)>;


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
            NearestNeighborLoop::Exterior { .. } 
            | NearestNeighborLoop::JointExterior { .. } => {
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
            _ => {
                for r in spans {
                    let p5 = *r.start() as NAIDX;
                    let p5s = p5 + 1;
                    let p3 = *r.end() as NAIDX;
                    let p3e = p3;
                    self.shift_iter(ltab, (p5, p5s), (p3, p3e), &loopdict, &mut neighbors);
                }
            },
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
                if (p == p5 && up + ltab.min_hairpin_size() >= uk) || 
                    (p == p3 && uk + ltab.min_hairpin_size() >= up) {
                        continue;
                } 
                if !ltab.can_pair(uk, up) {
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
 
    fn get_activation_energy<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> i32 {
        let (i_idx, o_idx, s_inner, s_outer) = self.get_loops(ltab, mv);

        let (_, inner_en) = ltab.get(i_idx);
        let (_, outer_en) = ltab.get(o_idx);

        let s_inner_en = ltab.energy_of_loop(&s_inner);
        let s_outer_en = ltab.energy_of_loop(&s_outer);

        (s_inner_en - inner_en)
            .max(s_outer_en - outer_en)
            .max((s_inner_en + s_outer_en) - (inner_en + outer_en))
    }

    fn get_loops<E: EnergyModel>(
        &self,
        ltab: &LoopTable<E>,
        mv: &Move,
    ) -> (usize, // init
          usize, // merge
          NearestNeighborLoop, // s_inner
          NearestNeighborLoop) // s_outer
    {
        match mv {
            Move::ShiftIK { i, j, k } => {
                let o_idx = ltab.loop_lookup(*i as usize);
                let i_idx = ltab.loop_lookup(*j as usize);
                let (outer, _en) = ltab.get(o_idx);
                let (inner, _en) = ltab.get(i_idx); 
                let combo = outer.join_loop(inner);
                let (p, q) = if *i < *k { (*i, *k) } else { (*k, *i) };
                let (s_outer, s_inner) = combo.split_loop(p, q);
                (i_idx, o_idx, s_outer, s_inner)
            },
            Move::ShiftJK { i, j, k } => {
                let o_idx = ltab.loop_lookup(*i as usize);
                let i_idx = ltab.loop_lookup(*j as usize);
                let (outer, _en) = ltab.get(o_idx);
                let (inner, _en) = ltab.get(i_idx); 
                let combo = outer.join_loop(inner);
                let (p, q) = if *j < *k { (*j, *k) } else { (*k, *j) };
                let (s_outer, s_inner) = combo.split_loop(p, q);
                (i_idx, o_idx, s_outer, s_inner)
            },
            _ => panic!("wrong move"),
        }

    }

    pub fn remove(
        &mut self, 
        index: &usize, 
    ) -> Option<Moves> {
        self.map
            .remove(index)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;


    #[test]
    fn test_three_way_setup() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GUACACGUCG");
        let pairings =       PairTable::try_from("..........").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut twn = ThreeWayNeighbors::default();
        let res = twn.compute_neighbors(&ltab, 0);
        assert!(
            res.is_empty(),
            "Expected no three-way shift neighbors."
        );
    }
 
    #[test]
    fn test_three_way_simple() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GUACCCCUCG");
        let pairings =       PairTable::try_from("...(.....)").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        println!("LT: {:?}", ltab);
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
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("CCCCUCG");
        let pairings =       PairTable::try_from("(.....)").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        println!("LT: {:?}", ltab);
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
