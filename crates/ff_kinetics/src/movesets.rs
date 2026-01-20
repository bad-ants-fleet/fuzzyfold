use std::fmt;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_energy::LoopDecomposition;
use ff_energy::NearestNeighborLoop;
use crate::Move;

type Moves = Vec<(Move, i32)>;
type IJMoves = Vec<(NAIDX, NAIDX, i32)>;
type LoopEntry = (NearestNeighborLoop, i32);

/// Stores all information about the current loop decomposition.
pub struct LoopTable<'a, E: EnergyModel> {
    sequence: &'a [Base],
    model: &'a E,
    loops: Vec<LoopEntry>,
    stale: Vec<usize>,
    loop_lookup: Vec<usize>,
    pair_lookup: IntMap<NAIDX, NAIDX>,
    energy: i32,
}

impl<'a, E: EnergyModel> LoopTable<'a, E> {
    pub fn energy(&self) -> i32 {
        self.energy
    }

    pub fn pair_lookup(&self, idx: &NAIDX) -> NAIDX {
        self.pair_lookup[idx]
    }

    pub fn geti(&self, i: usize) -> &LoopEntry {
        &self.loops[self.loop_lookup[i]]
    }

    pub fn get(&self, idx: usize) -> &LoopEntry {
        self.loops
            .get(idx)
            .expect("Invalid LoopCache index")
    }

    pub fn insert_loop(&mut self, nn_loop: NearestNeighborLoop, nn_energy: i32) -> usize {
        if let Some(idx) = self.stale.pop() {
            self.loops[idx] = (nn_loop, nn_energy);
            idx
        } else {
            let idx = self.loops.len();
            self.loops.push((nn_loop, nn_energy));
            idx
        }
    }
}

impl<'a, E: EnergyModel> fmt::Display for LoopTable<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Convert sequence to string
        let mut dbr = vec!['.'; self.loop_lookup.len()];
        for (i, j) in &self.pair_lookup {
            dbr[*i as usize] = '(';
            dbr[*j as usize] = ')';
        }
        let dbr_str: String = dbr.into_iter().collect();
        write!(f, "{}", dbr_str)
    }
}

impl<'a, E: EnergyModel> fmt::Debug for LoopTable<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopTable")
            .field("loop_lookup", &self.loop_lookup)
            .field("pair_lookup", &self.pair_lookup)
            .finish()
    }
}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a [Base], &T, &'a E)> 
for LoopTable<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a E)
    ) -> Result<Self, Self::Error> {

        let mut loops = Vec::new();
        let mut loop_lookup: Vec<usize> = vec![0; sequence.len()];
        let mut pair_lookup: IntMap<NAIDX, NAIDX>  = IntMap::default();
        let mut energy = 0;

        pairings.for_each_loop(|l| {
            let loop_energy = model.energy_of_loop(sequence, l);
            energy += loop_energy;

            if let Some((i, j)) = l.closing() {
                pair_lookup.insert(i as NAIDX, j as NAIDX); 
            }

            let loop_index = loops.len();
            for k in l.inclusive_unpaired_indices() {
                loop_lookup[k] = loop_index;
            }
            loops.push((l.to_owned(), loop_energy));
        });

        Ok(LoopTable {
            sequence,
            model,
            loops,
            stale: Vec::new(),
            loop_lookup,
            pair_lookup,
            energy,
        })
    }
}

pub struct AddDelMoves<'a, E: EnergyModel> {
    loop_table: LoopTable<'a, E>,
    /// registry index to list of (i, j, deltaE)
    add_neighbors: IntMap<usize, IJMoves>,
    /// pair id to deltaE
    del_neighbors: IntMap<NAIDX, i32>, 
}

impl<'a, E: EnergyModel> fmt::Display for AddDelMoves<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loop_table)
    }
}

impl<'a, E: EnergyModel> AddDelMoves<'a, E> {
    pub fn loop_table(&self) -> &LoopTable<'a, E> {
        &self.loop_table
    }

    pub fn add_neighbors(&self) -> &IntMap<usize, IJMoves> {
        &self.add_neighbors
    }

    pub fn del_neighbors(&self) -> &IntMap<NAIDX, i32> {
        &self.del_neighbors
    }

    fn init_del_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for (&i, &j) in ltab.pair_lookup.iter() {
            let delta = self.eval_del_move(i, j);
            self.del_neighbors.insert(i, delta);
        }
    }

    fn update_del_moves(
        &mut self,
        pairs: &[(NAIDX, NAIDX)],
        cap: usize,
    ) -> Moves
    {
        let mut change = Vec::with_capacity(pairs.len() + cap);
        for &(i, j) in pairs {
            let delta = self.eval_del_move(i, j);
            self.del_neighbors.insert(i, delta);
            change.push((Move::Del{ i, j }, delta));
        }
        change
    }

    /// Returns how the free energy changes if the move is applied.
    fn eval_del_move(&self, i: NAIDX, j: NAIDX) -> i32 {
        let ltab = &self.loop_table;
        let (outer, o_en) = ltab.geti(i as usize);
        let (inner, i_en) = ltab.geti(j as usize);
        let combo = outer.join_loop(inner);
        let combo_energy = ltab.model.energy_of_loop(ltab.sequence, &combo);
        combo_energy - (o_en + i_en)
    }

    fn init_add_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for lli in 0..ltab.loops.len() {
            let neighbors = self.get_add_neighbors_per_loop(lli);
            self.add_neighbors.insert(lli, neighbors);
        }
    }

    fn get_add_neighbors_per_loop(&self, index: usize) -> IJMoves {
        let ltab = &self.loop_table;
        let (combo, energy) = ltab.get(index);
        let unpaired = combo.unpaired_indices();

        let mut neighbors = Vec::new(); 
        for (idx_i, &i) in unpaired.iter().enumerate() {
            for &j in &unpaired[idx_i + 1..] {
                if j <= i + ltab.model.min_hairpin_size() {
                    continue;
                }
                if ltab.model.can_pair(ltab.sequence[i], ltab.sequence[j]) {
                    let i = i as NAIDX;
                    let j = j as NAIDX;
                    let (outer, inner) = combo.split_loop(i, j);
                    let outer_energy = ltab.model.energy_of_loop(ltab.sequence, &outer);
                    let inner_energy = ltab.model.energy_of_loop(ltab.sequence, &inner);
                    // How does the free energy change if the move is applied.
                    let delta = (outer_energy + inner_energy) - energy;
                    neighbors.push((i, j, delta)); 
                }
            }
        }
        neighbors
    }
 
    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let delta = self.del_neighbors.remove(&i).expect("Missing pair_list entry.");

        let ltab = &mut self.loop_table;
        if j != ltab.pair_lookup.remove(&i).expect("Missing pair_list entry.") {
            panic!("Inconsistent pair-list entry.");
        }
        let o_index = ltab.loop_lookup[i as usize];
        let i_index = ltab.loop_lookup[j as usize];

        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);
        let combo = outer.join_loop(inner);
        let combo_energy = (o_en + i_en) + delta;
        debug_assert_eq!(combo_energy, 
            ltab.model.energy_of_loop(ltab.sequence, &combo));

        for k in combo.inclusive_unpaired_indices() {
            debug_assert!(ltab.loop_lookup[k] == o_index || ltab.loop_lookup[k] == i_index);
            ltab.loop_lookup[k] = o_index;
        }

        let cpairs = &combo.pairs();

        ltab.loops[o_index] = (combo, combo_energy);
        ltab.stale.push(i_index);
        ltab.energy += delta;
 
        // Those include the neighbors
        let update_add_moves = self.get_add_neighbors_per_loop(o_index);
        let mut new_moves = self.update_del_moves(cpairs, update_add_moves.len());
        new_moves.extend(
            update_add_moves.iter()
            .map(|&(i, j, d)| (Move::Add { i, j }, d)));

        self.add_neighbors.insert(o_index, update_add_moves);
        self.add_neighbors.remove(&i_index).expect("at least empty list.");

        (Vec::new(), new_moves)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;

        let c_index = ltab.loop_lookup[i as usize];
        debug_assert_eq!(c_index, ltab.loop_lookup[j as usize], "Missing loop_lookup entry for j.");

        let (combo, c_en) = ltab.get(c_index).clone();
        let old_add_moves = self.add_neighbors.remove(&c_index).expect("Old combo moves")
            .into_iter()
            .filter(|&(p, q, _)| {
                !(q < i || j < p || (i < p && q < j) || (p < i && j < q))
            })
            .map(|(i, j, d)| (Move::Add { i, j }, d))
            .collect();

        let (outer, inner) = combo.split_loop(i, j);
        
        let o_en = ltab.model.energy_of_loop(ltab.sequence, &outer);
        ltab.loops[c_index] = (outer, o_en);
        let i_en = ltab.model.energy_of_loop(ltab.sequence, &inner);
        let i_index = ltab.insert_loop(inner, i_en);

        let delta = (o_en + i_en) - c_en;
        ltab.energy += delta;
        ltab.pair_lookup.insert(i, j);
        let (outer, _) = ltab.get(c_index);
        for k in outer.inclusive_unpaired_indices() {
            ltab.loop_lookup[k] = c_index;
        }
        let (inner, _) = ltab.get(i_index);
        for k in inner.inclusive_unpaired_indices() {
            ltab.loop_lookup[k] = i_index;
        }
        let outer_add_moves = self.get_add_neighbors_per_loop(c_index);
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index);

        let cap = outer_add_moves.len() + inner_add_moves.len() + 1;
        let mut new_moves = self.update_del_moves(&combo.pairs(), cap);
        new_moves.push((Move::Del{ i, j }, -delta));
        new_moves.extend(
            outer_add_moves.iter()
            .chain(inner_add_moves.iter())
            .map(|&(i, j, d)| (Move::Add { i, j }, d)));

        self.add_neighbors.insert(c_index, outer_add_moves);
        self.add_neighbors.insert(i_index, inner_add_moves);
        self.del_neighbors.insert(i, -delta);
        (old_add_moves, new_moves)
    }

}

impl<'a, E: EnergyModel> From<LoopTable<'a, E>> for AddDelMoves<'a, E> {
    fn from(loop_table: LoopTable<'a, E>) -> Self {
        let mut moves = AddDelMoves {
            loop_table,
            add_neighbors: IntMap::default(),
            del_neighbors: IntMap::default(),
        };
        moves.init_add_neighbors();
        moves.init_del_neighbors();
        moves
    }
}

impl<'a, E: EnergyModel> From<&LoopTable<'a, E>> for DotBracketVec {
    fn from(ltab: &LoopTable<'a, E>) -> Self {
        // Use the same logic as your Display impl, but avoid allocating a String unnecessarily
        let mut vec = vec![DotBracket::Unpaired; ltab.loop_lookup.len()];
        for (i, j) in &ltab.pair_lookup {
            vec[*i as usize] = DotBracket::Open;
            vec[*j as usize] = DotBracket::Close;
        }
        DotBracketVec(vec)
    }
}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a [Base], &T, &'a E)> 
for AddDelMoves<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a E)
    ) -> Result<Self, Self::Error> {

        let mut loops = Vec::new();
        let mut loop_lookup: Vec<usize> = vec![0; sequence.len()];
        let mut pair_lookup: IntMap<NAIDX, NAIDX>  = IntMap::default();
        let mut energy = 0;

        pairings.for_each_loop(|l| {
            let loop_energy = model.energy_of_loop(sequence, l);
            energy += loop_energy;

            if let Some((i, j)) = l.closing() {
                pair_lookup.insert(i as NAIDX, j as NAIDX); 
            }

            let loop_index = loops.len();
            for k in l.inclusive_unpaired_indices() {
                loop_lookup[k] = loop_index;
            }
            loops.push((l.to_owned(), loop_energy));
        });

        Ok(AddDelMoves::from(LoopTable {
            sequence,
            model,
            loops,
            stale: Vec::new(),
            loop_lookup,
            pair_lookup,
            energy,
        }))
    }
}



//impl<'a, E: EnergyModel> LoopCache<'a, E> {
//
//    fn get_shift_loop_neighbors(&self, index: usize, ll: &[usize]) -> Vec<(Move, i32)> {
//        let (combo, combo_energy) = self.get(index);
//        let mut neighbors = Vec::new(); 
//
//        match combo {
//            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
//                let ui = *i as usize;
//                let uj = *j as usize;
//                let (outer, outer_energy) = self.get(ll[ui]); 
//                let shift_combo = outer.join_loop(combo);
//                let energy = combo_energy + outer_energy;
//
//                for k in ui+1..uj {
//                    if k > ui + self.model.min_hairpin_size() && 
//                        self.model.can_pair(self.sequence[ui], self.sequence[k]) {
//                            let (s_outer, s_inner) = shift_combo.split_loop(*i, k as NAIDX);
//                            let s_outer_energy = self.model.energy_of_loop(self.sequence, &s_outer);
//                            let s_inner_energy = self.model.energy_of_loop(self.sequence, &s_inner);
//                            let delta = (s_outer_energy + s_inner_energy) - energy;
//                            neighbors.push((Move::ShiftI { i: *i, j: *j, k: k as NAIDX }, delta));
//                    }
//                    if k < uj - self.model.min_hairpin_size() && 
//                        self.model.can_pair(self.sequence[k], self.sequence[uj]) {
//                            let (s_outer, s_inner) = shift_combo.split_loop(k as NAIDX, *j);
//                            let s_outer_energy = self.model.energy_of_loop(self.sequence, &s_outer);
//                            let s_inner_energy = self.model.energy_of_loop(self.sequence, &s_inner);
//                            let delta = (s_outer_energy + s_inner_energy) - energy;
//                            neighbors.push((Move::ShiftJ { i: *i, j: *j, k: k as NAIDX }, delta));
//                    }
//                }
//            }
//            _ => {warn!("not implemented")},
//        }
//        neighbors
//    }
//}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Walker;

    #[test]
    fn test_add_then_del_roundtrip() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GCCCCGGUCA");
        let pairings =       PairTable::try_from("..........").unwrap();

        let ltab = LoopTable::try_from((&sequence[..], &pairings, &model)).unwrap();
        let mut adm = AddDelMoves::from(ltab);

        // Clone neighbor list so we don’t mutate while iterating
        let neighbors: Vec<(Move, i32)> = adm.propose_moves().collect();

        for (mv, de) in neighbors {
            let en0 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv, de, en0);

            // add pair
            let _ = adm.apply_move(&mv);
            let en1 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv, de, en1);

            let _ = adm.apply_move(&mv.inverse());
            let en2 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv.inverse(), de, en2);

            assert_eq!(en0, en2, "roundtrip energy mismatch");
        }
    }

    #[test]
    fn test_add_then_del_bug() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GCCCCGGUCA");
        let pairings1 =      PairTable::try_from("((....).).").unwrap();
        let pairings2 =      PairTable::try_from("..........").unwrap();

        let ltab1 = LoopTable::try_from((&sequence[..], &pairings1, &model)).unwrap();
        let ltab2 = LoopTable::try_from((&sequence[..], &pairings2, &model)).unwrap();

        let adm1 = AddDelMoves::from(ltab1);
        let mut adm2 = AddDelMoves::from(ltab2);

        let nb1: Vec<_> = adm1.propose_moves().collect();
        println!("{}: {:?}", adm1.current_energy(), nb1);

        let _ = adm2.apply_add_move(0, 8);
        let _ = adm2.apply_add_move(1, 6);
        let nb2: Vec<_> = adm2.propose_moves().collect();
        println!("{}: {:?}", adm2.current_energy(), nb2);
        assert_eq!(nb1, nb2);
    }

}

