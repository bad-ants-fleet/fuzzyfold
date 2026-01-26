use std::fmt;
use ff_energy::NucleotideVec;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::EnergyModel;
use ff_energy::LoopDecomposition;
use crate::Move;
use crate::movesets::LoopTable;

type Moves = Vec<(Move, i32)>;

pub struct AddDelMoves<'a, E: EnergyModel> {
    loop_table: LoopTable<'a, E>,
    /// registry index to list of (i, j, deltaE)
    add_neighbors: IntMap<usize, Moves>,
    /// pair id to deltaE
    del_neighbors: IntMap<NAIDX, i32>, 
}

impl<'a, E: EnergyModel> Clone for AddDelMoves<'a, E> {
    fn clone(&self) -> Self {
        Self {
            loop_table: self.loop_table.clone(),   // clones refs + vectors
            add_neighbors: self.add_neighbors.clone(),
            del_neighbors: self.del_neighbors.clone(),
        }
    }
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

    pub fn add_neighbors(&self) -> &IntMap<usize, Moves> {
        &self.add_neighbors
    }

    pub fn del_neighbors(&self) -> &IntMap<NAIDX, i32> {
        &self.del_neighbors
    }

    fn init_del_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for (&i, &j) in ltab.pairs() {
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
        let combo_energy = ltab.energy_of_loop(&combo);
        combo_energy - (o_en + i_en)
    }

    fn init_add_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for lli in 0..ltab.loops_len() {
            let neighbors = self.get_add_neighbors_per_loop(lli);
            self.add_neighbors.insert(lli, neighbors);
        }
    }

    fn get_add_neighbors_per_loop(&self, index: usize) -> Moves {
        let ltab = &self.loop_table;
        let (combo, energy) = ltab.get(index);
        let unpaired = combo.unpaired_indices();

        let mut neighbors = Vec::new(); 
        for (idx_i, &i) in unpaired.iter().enumerate() {
            for &j in &unpaired[idx_i + 1..] {
                if j <= i + ltab.min_hairpin_size() {
                    continue;
                }
                if ltab.can_pair(i, j) {
                    let i = i as NAIDX;
                    let j = j as NAIDX;
                    let (outer, inner) = combo.split_loop(i, j);
                    let outer_energy = ltab.energy_of_loop(&outer);
                    let inner_energy = ltab.energy_of_loop(&inner);
                    // How does the free energy change if the move is applied.
                    let delta = (outer_energy + inner_energy) - energy;
                    neighbors.push((Move::Add { i, j }, delta)); 
                }
            }
        }
        neighbors
    }
 
    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;

        let o_index = ltab.loop_lookup(i as usize);
        let i_index = ltab.loop_lookup(j as usize);
        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);

        let delta = self.del_neighbors.remove(&i)
            .expect("Missing del neighbor.");
        let _ = self.add_neighbors.remove(&i_index)
            .expect("Old add moves");

        let combo = outer.join_loop(inner);
        let combo_pairs = &combo.pairs();
        let c_en = (o_en + i_en) + delta;
        debug_assert_eq!(c_en, ltab.energy_of_loop(&combo));

        let c_index = ltab.insert_loopentry(Some(o_index), (combo, c_en));
        ltab.mark_stale(i_index);
        ltab.update_lookup(c_index);
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }

        // Those include the neighbors
        let update_add_moves = self.get_add_neighbors_per_loop(c_index);
        let mut new_moves = self.update_del_moves(combo_pairs, update_add_moves.len());
        new_moves.extend(update_add_moves.iter().cloned());

        self.add_neighbors.insert(o_index, update_add_moves);

        (vec![(Move::Del { i, j }, delta)], new_moves)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;

        let c_index = ltab.loop_lookup(i as usize);
        debug_assert_eq!(c_index, ltab.loop_lookup(j as usize), "Missing loop_lookup entry for j.");
        let (combo, c_en) = ltab.get(c_index);
        let combo_pairs = &combo.pairs();
        let old_add_moves = self.add_neighbors
            .remove(&c_index)
            .expect("Old combo moves")
            .into_iter()
            .filter(|(amv, _)| amv.conflicts((i, j)))
            .collect();

        let (outer, inner) = combo.split_loop(i, j);
        let o_en = ltab.energy_of_loop(&outer);
        let i_en = ltab.energy_of_loop(&inner);
        let delta = (o_en + i_en) - c_en;

        ltab.insert_pair(i, j);
        let o_index = ltab.insert_loopentry(Some(c_index), (outer, o_en));
        let i_index = ltab.insert_loopentry(None, (inner, i_en));
        ltab.update_lookup(o_index);
        ltab.update_lookup(i_index);


        let outer_add_moves = self.get_add_neighbors_per_loop(o_index);
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index);
        let cap = outer_add_moves.len() + inner_add_moves.len() + 1;
        let mut new_moves = self.update_del_moves(combo_pairs, cap);
        new_moves.push((Move::Del{ i, j }, -delta));
        new_moves.extend([
            &outer_add_moves,
            &inner_add_moves,
        ].into_iter()
        .flat_map(|v| v.iter().cloned()));

        self.add_neighbors.insert(o_index, outer_add_moves);
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

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a NucleotideVec, &T, &'a E)> 
for AddDelMoves<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a NucleotideVec, &T, &'a E)
    ) -> Result<Self, Self::Error> {
        Ok(AddDelMoves::from(LoopTable::try_from((sequence, pairings, model))?))
    }
}

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

        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
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

        let ltab1 = LoopTable::try_from((&sequence, &pairings1, &model)).unwrap();
        let ltab2 = LoopTable::try_from((&sequence, &pairings2, &model)).unwrap();

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

    #[test]
    fn test_development_bug02() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GUACACGUCG");
        let pairings =       PairTable::try_from("..........").unwrap();

        let mut adm = AddDelMoves::try_from((&sequence, &pairings, &model)).unwrap();
        let en1 = adm.current_energy();
        let nb1: Vec<_> = adm.propose_moves().collect();
        for (mv, d) in nb1 {
            println!("pp {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::Add { i: 0, j: 8 });
        println!("Applied add");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::Del { i: 0, j: 8 });
        println!("Applied del");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }
        assert_eq!(en1, adm.current_energy());
    }
}

