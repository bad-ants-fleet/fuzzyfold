use std::fmt;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_energy::LoopDecomposition;

use crate::Move;
use crate::movesets::LoopTable;
use crate::movesets::three_way_shifts::ThreeWayNeighbors;
use crate::movesets::four_way_shifts::FourWayNeighbors;

type Moves = Vec<(Move, i32)>;

#[derive(Debug, Clone, Copy)]
pub struct ShiftConfig {
    pub three_way: bool,
    pub four_way: bool,
}

pub struct AddDelShiftMoves<'a, E: EnergyModel> {
    loop_table: LoopTable<'a, E>,
    add_neighbors: IntMap<usize, Moves>,
    del_neighbors: IntMap<NAIDX, i32>, 
    three_way_shift_neighbors: Option<ThreeWayNeighbors>,
    four_way_shift_neighbors: Option<FourWayNeighbors>,
}

impl<'a, E: EnergyModel> AddDelShiftMoves<'a, E> {
    pub fn loop_table(&self) -> &LoopTable<'a, E> {
        &self.loop_table
    }

    pub fn add_neighbors(&self) -> &IntMap<usize, Moves> {
        &self.add_neighbors
    }

    pub fn del_neighbors(&self) -> &IntMap<NAIDX, i32> {
        &self.del_neighbors
    }

    pub fn three_way_shift_neighbors(&self) -> &Option<ThreeWayNeighbors> {
        &self.three_way_shift_neighbors
    }

    pub fn four_way_shift_neighbors(&self) -> &Option<FourWayNeighbors> {
        &self.four_way_shift_neighbors
    }

    /// Activation energy -> for add moves it is delta-E
    /// We should be able to change this for alternative 
    /// rate models, but do not break detailed balance.
    fn get_add_activation_energy(&self, 
        index: usize,
        i: NAIDX, 
        j: NAIDX,
    ) -> i32 {
        let ltab = &self.loop_table;
        let (combo, combo_energy) = ltab.get(index);
        let (outer, inner) = combo.split_loop(i, j);
        let outer_energy = ltab.energy_of_loop(&outer);
        let inner_energy = ltab.energy_of_loop(&inner);
        (outer_energy + inner_energy) - combo_energy
    }

    /// Returns how the free energy changes if the move is applied.
    fn get_del_activation_energy(&self, 
        i: NAIDX, 
        j: NAIDX,
    ) -> i32 {
        let ltab = &self.loop_table;
        let (outer, o_en) = ltab.geti(i as usize);
        let (inner, i_en) = ltab.geti(j as usize);
        let combo = outer.join_loop(inner);
        let combo_energy = ltab.energy_of_loop(&combo);
        combo_energy - (o_en + i_en)
    }

    fn init_del_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for (&i, &j) in ltab.pairs() {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
        }
    }

    fn init_loop_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for lli in 0..ltab.loops_len() {
            self.get_add_neighbors_per_loop(lli);
            let ltab = &self.loop_table;
            if let Some(tw) = &mut self.three_way_shift_neighbors {
                tw.compute_neighbors(ltab, lli);
            }
            if let Some(fw) = &mut self.four_way_shift_neighbors {
                fw.compute_neighbors(ltab, lli);
            }
        }
    }

    fn get_add_neighbors_per_loop(&mut self, index: usize) -> &Moves {
        let ltab = &self.loop_table;
        let (combo, _) = ltab.get(index);
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
                    let barrier = self.get_add_activation_energy(index, i, j);
                    neighbors.push((Move::Add { i, j }, barrier)); 
                }
            }
        }
        self.add_neighbors.insert(index, neighbors);
        &self.add_neighbors[&index]
    }

    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;
        let o_index = ltab.loop_lookup(i as usize);
        let i_index = ltab.loop_lookup(j as usize);
        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);

        let delta = self.del_neighbors
            .remove(&i)
            .expect("Missing old del move.");
        // those are deleted because we don't reuse the index later!
        // All of these moves will be valid in the future, but their 
        // energy evaluation changes.
        let _ = self.add_neighbors
            .remove(&i_index)
            .expect("Missing old add moves");

        let old_tw_shift_outer = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&o_index));
        let old_tw_shift_inner = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&i_index));
        let old_fw_shift_outer = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&o_index));
        let old_fw_shift_inner = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&i_index));

        let old_moves: Moves =
            [
            old_tw_shift_outer,
            old_tw_shift_inner,
            old_fw_shift_outer,
            old_fw_shift_inner,
            ]
            .into_iter()
            .flatten()      // Option<Moves> → Moves
            .flatten()      // Moves → (Move, i32)
            .filter(move |(mv, _)| mv.conflicts((i, j)))
            .chain(std::iter::once((Move::Del { i, j }, delta)))
            .collect();

        let combo = outer.join_loop(inner);
        let combo_pairs = combo.pairs().to_vec();
        let c_en = o_en + i_en + delta;
        debug_assert_eq!(c_en, ltab.energy_of_loop(&combo));

        // Update the loop table with all new data.
        let c_index = ltab.insert_loopentry(Some(o_index), (combo, c_en));
        ltab.mark_stale(i_index);
        ltab.update_lookup(c_index);
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }

        // Those include the neighbors
        let new_add_moves: Moves = self.get_add_neighbors_per_loop(c_index).clone();

        let ltab = &self.loop_table;
        let new_tw_shift_moves: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, c_index).clone())
            .unwrap_or_default();
        let new_fw_shift_moves: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, c_index).clone())
            .unwrap_or_default();
        let cap = new_add_moves.len() 
            + new_tw_shift_moves.len()
            + new_fw_shift_moves.len()
            + combo_pairs.len() + 1;
        let mut new_moves = Vec::with_capacity(cap);
        for (i, j) in combo_pairs {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        new_moves.extend(new_add_moves);
        new_moves.extend(new_tw_shift_moves);
        new_moves.extend(new_fw_shift_moves);

        (old_moves, new_moves)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;

        // Get the original "combo" loop & the conflicting moves.
        let c_index = ltab.loop_lookup(i as usize);
        debug_assert_eq!(c_index, ltab.loop_lookup(j as usize), "Missing loop_lookup entry for j.");
        let (combo, c_en) = ltab.get(c_index);
        let combo_pairs = &combo.pairs();


        let old_add_moves = self.add_neighbors.remove(&c_index);
        let old_tw_shift = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&c_index));
        let old_fw_shift = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&c_index));

        let old_moves: Moves =
            [
            old_add_moves,
            old_tw_shift,
            old_fw_shift,
            ]
            .into_iter()
            .flatten()      // Option<Moves> → Moves
            .flatten()      // Moves → (Move, i32)
            .filter(move |(mv, _)| mv.conflicts((i, j)))
            .collect();

        // Calculate the energies for the new loops (again...).
        let (outer, inner) = combo.split_loop(i, j);
        let o_en = ltab.energy_of_loop(&outer);
        let i_en = ltab.energy_of_loop(&inner);
        let delta = (o_en + i_en) - c_en;

        // Update the loop table with all new data.
        ltab.insert_pair(i, j);
        let o_index = ltab.insert_loopentry(Some(c_index), (outer, o_en));
        let i_index = ltab.insert_loopentry(None, (inner, i_en));
        ltab.update_lookup(o_index);
        ltab.update_lookup(i_index);

        let outer_add_moves = self.get_add_neighbors_per_loop(o_index).clone();
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index).clone();

        let ltab = &self.loop_table;
        let tw_shift_moves_outer: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o_index).clone())
            .unwrap_or_default();
        let tw_shift_moves_inner: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, i_index).clone())
            .unwrap_or_default();
        let fw_shift_moves_outer: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o_index).clone())
            .unwrap_or_default();
        let fw_shift_moves_inner: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, i_index).clone())
            .unwrap_or_default();
        let cap = outer_add_moves.len() 
            + inner_add_moves.len() 
            + outer_add_moves.len()
            + tw_shift_moves_inner.len()
            + tw_shift_moves_outer.len() 
            + fw_shift_moves_inner.len()
            + fw_shift_moves_outer.len() 
            + combo_pairs.len() + 1;

        let mut new_moves = Vec::with_capacity(cap);
        for &(i, j) in combo_pairs {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        self.del_neighbors.insert(i, -delta);
        new_moves.push((Move::Del{ i, j }, -delta));
        new_moves.extend([
            &outer_add_moves,
            &inner_add_moves,
            &tw_shift_moves_outer,
            &tw_shift_moves_inner,
            &fw_shift_moves_outer,
            &fw_shift_moves_inner,
        ].into_iter() 
        .flat_map(|v| v.iter().cloned()));

        (old_moves, new_moves)
    }

    pub fn apply_three_way_shift_move(
        &mut self, 
        mv: &Move, 
        k: NAIDX,
    ) -> (Moves, Moves) {
        let ltab = &mut self.loop_table;
        let (i, j) = mv.deleted_pair();
        let (p, q) = mv.added_pair();

        let o_index = ltab.loop_lookup(i as usize);
        let i_index = ltab.loop_lookup(j as usize);
        let k_index = ltab.loop_lookup(k as usize);
        debug_assert!(o_index == k_index || i_index == k_index);

        let delta = self.del_neighbors.remove(&i)
            .expect("Missing del neighbor.");

        let old_add_init = self.add_neighbors.remove(&k_index);
        let old_tw_shift_outer = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&o_index));
        let old_tw_shift_inner = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&i_index));
        let old_fw_shift_outer = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&o_index));
        let old_fw_shift_inner = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&i_index));

        let old_moves: Moves = old_add_init.into_iter()
            .flatten()
            .filter(|(mv, _)| mv.conflicts((p, q)))
            .chain(old_tw_shift_outer.into_iter().flatten()) 
            .chain(old_tw_shift_inner.into_iter().flatten())
            .chain(old_fw_shift_outer.into_iter().flatten())
            .chain(old_fw_shift_inner.into_iter().flatten())
            .filter(|(mv, _)| mv.conflicts((p, q)) || mv.conflicts((i, j)))
            .chain(std::iter::once((Move::Del { i, j }, delta)))
            .collect();

        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);
        let combo = outer.join_loop(inner);
        let combo_pairs = &combo.pairs();
        let c_en = o_en + i_en + delta;
        debug_assert_eq!(c_en, ltab.energy_of_loop(&combo));
        let (new_outer, new_inner) = combo.split_loop(p, q);
        let new_o_en = ltab.energy_of_loop(&new_outer);
        let new_i_en = ltab.energy_of_loop(&new_inner);
        let del_delta = c_en - (new_o_en + new_i_en);

        let _ = ltab.insert_loopentry(Some(o_index), (new_outer, new_o_en));
        let _ = ltab.insert_loopentry(Some(i_index), (new_inner, new_i_en));
        ltab.update_lookup(o_index);
        ltab.update_lookup(i_index);
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        ltab.insert_pair(p, q);

        let outer_add_moves = self.get_add_neighbors_per_loop(o_index).clone();
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index).clone();

        let ltab = &self.loop_table;
        let tw_shift_moves_outer: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o_index).clone())
            .unwrap_or_default();
        let tw_shift_moves_inner: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, i_index).clone())
            .unwrap_or_default();
        let fw_shift_moves_outer: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o_index).clone())
            .unwrap_or_default();
        let fw_shift_moves_inner: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, i_index).clone())
            .unwrap_or_default();
        let cap = outer_add_moves.len() 
            + inner_add_moves.len() 
            + outer_add_moves.len()
            + tw_shift_moves_inner.len()
            + tw_shift_moves_outer.len() 
            + fw_shift_moves_inner.len()
            + fw_shift_moves_outer.len() 
            + combo_pairs.len() + 1;

        let mut new_moves = Vec::with_capacity(cap);
        for &(i, j) in combo_pairs {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        new_moves.extend([
            &outer_add_moves,
            &inner_add_moves,
            &tw_shift_moves_outer,
            &tw_shift_moves_inner,
            &fw_shift_moves_outer,
            &fw_shift_moves_inner,
        ].into_iter() 
        .flat_map(|v| v.iter().cloned()));
        new_moves.push((Move::Del { i: p, j: q }, del_delta));

        self.del_neighbors.insert(p, del_delta);

        (old_moves, new_moves)
    }

    pub fn apply_four_way_shift_move(
        &mut self, 
        mv: &Move, 
    ) -> (Moves, Moves) {
        let ltab = &mut self.loop_table;

        let ((i, j), (k, l)) = mv.deleted_pairs();
        let deld1 = self.del_neighbors.remove(&i).expect("Missing del neighbor.");
        let deld2 = self.del_neighbors.remove(&k).expect("Missing del neighbor.");

        let (it_idx, m1_idx, m2_idx, inner0, outer1, outer2) = 
            self.four_way_shift_neighbors.as_ref().expect("applying four way!")
            .get_loops(ltab, mv);

        let old_add_init = self.add_neighbors.remove(&it_idx);
        let old_add_merge1 = self.add_neighbors.remove(&m1_idx);
        let old_add_merge2 = self.add_neighbors.remove(&m2_idx);
        let old_tw_shift_init = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&it_idx));
        let old_tw_shift_merge1 = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&m1_idx));
        let old_tw_shift_merge2 = self.three_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&m2_idx));
        let old_fw_shift_init = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&it_idx));
        let old_fw_shift_merge1 = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&m1_idx));
        let old_fw_shift_merge2 = self.four_way_shift_neighbors.as_mut()
            .and_then(|map| map.remove(&m2_idx));

        let ((p, q), (m, n)) = mv.added_pairs();
        let old_moves: Moves = old_add_init.into_iter()
            .chain(old_add_merge1) 
            .chain(old_add_merge2) 
            .chain(old_tw_shift_init) 
            .chain(old_tw_shift_merge1) 
            .chain(old_tw_shift_merge2) 
            .chain(old_fw_shift_init) 
            .chain(old_fw_shift_merge1) 
            .chain(old_fw_shift_merge2) 
            .flatten()
            .filter(|(mv, _)| mv.conflicts((p, q)) || mv.conflicts((m, n)))
            .chain(std::iter::once((Move::Del { i, j }, deld1)))
            .chain(std::iter::once((Move::Del { i: k, j: l }, deld2)))
            .collect();

        let inner0_pairs = inner0.pairs().to_vec();
        let outer1_pairs = outer1.pairs().to_vec();
        let outer2_pairs = outer2.pairs().to_vec();
        let inner0_en = ltab.energy_of_loop(&inner0);  
        let outer1_en = ltab.energy_of_loop(&outer1);
        let outer2_en = ltab.energy_of_loop(&outer2);

        let in_idx = ltab.insert_loopentry(Some(it_idx), (inner0, inner0_en));
        let o1_idx = ltab.insert_loopentry(Some(m1_idx), (outer1, outer1_en));
        let o2_idx = ltab.insert_loopentry(Some(m2_idx), (outer2, outer2_en));
        ltab.update_lookup(in_idx);
        ltab.update_lookup(o1_idx);
        ltab.update_lookup(o2_idx);
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        if l != ltab.delete_pair(&k) {
            panic!("Inconsistent pair-list entry.");
        }
        ltab.insert_pair(p, q);
        ltab.insert_pair(m, n);

        let inner0_add_moves = self.get_add_neighbors_per_loop(in_idx).clone();
        let outer1_add_moves = self.get_add_neighbors_per_loop(o1_idx).clone();
        let outer2_add_moves = self.get_add_neighbors_per_loop(o2_idx).clone();


        let ltab = &self.loop_table;
        let tw_shift_moves_inner0: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, in_idx).clone())
            .unwrap_or_default();
        let tw_shift_moves_outer1: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o1_idx).clone())
            .unwrap_or_default();
        let tw_shift_moves_outer2: Moves = self
            .three_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o2_idx).clone())
            .unwrap_or_default();
        let fw_shift_moves_inner0: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, in_idx).clone())
            .unwrap_or_default();
        let fw_shift_moves_outer1: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o1_idx).clone())
            .unwrap_or_default();
        let fw_shift_moves_outer2: Moves = self
            .four_way_shift_neighbors.as_mut()
            .map(|map| map.compute_neighbors(ltab, o2_idx).clone())
            .unwrap_or_default();

       let cap = inner0_add_moves.len() 
            + outer1_add_moves.len()
            + outer2_add_moves.len()
            + tw_shift_moves_inner0.len()
            + tw_shift_moves_outer1.len()
            + tw_shift_moves_outer2.len()
            + fw_shift_moves_inner0.len()
            + fw_shift_moves_outer1.len()
            + fw_shift_moves_outer2.len()
            + inner0_pairs.len() 
            + outer1_pairs.len()
            + outer2_pairs.len();

        let mut new_moves = Vec::with_capacity(cap);
        for (i, j) in inner0_pairs {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        for &(i, j) in &outer1_pairs {
            if i == p || i == m {
                continue;
            }
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        for &(i, j) in &outer2_pairs {
            if i == p || i == m {
                continue;
            }
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            new_moves.push((Move::Del{ i, j }, delta));
        }
        new_moves.extend([
            &inner0_add_moves,
            &outer1_add_moves,
            &outer2_add_moves,
            &tw_shift_moves_inner0,
            &tw_shift_moves_outer1,
            &tw_shift_moves_outer2,
            &fw_shift_moves_inner0,
            &fw_shift_moves_outer1,
            &fw_shift_moves_outer2,
        ].into_iter() 
        .flat_map(|v| v.iter().cloned()));

        (old_moves, new_moves)
    }
}

impl<'a, E: EnergyModel> Clone for AddDelShiftMoves<'a, E> {
    fn clone(&self) -> Self {
        Self {
            loop_table: self.loop_table.clone(),   // clones refs + vectors
            add_neighbors: self.add_neighbors.clone(),
            del_neighbors: self.del_neighbors.clone(),
            three_way_shift_neighbors: self.three_way_shift_neighbors.clone(),
            four_way_shift_neighbors: self.four_way_shift_neighbors.clone(),
        }
    }
}

impl<'a, E: EnergyModel> fmt::Display for AddDelShiftMoves<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loop_table)
    }
}

impl<'a, E: EnergyModel> From<(LoopTable<'a, E>, ShiftConfig)>
    for AddDelShiftMoves<'a, E>
{
    fn from((loop_table, cfg): (LoopTable<'a, E>, ShiftConfig)) -> Self {
        let mut moves = AddDelShiftMoves {
            loop_table,
            add_neighbors: IntMap::default(),
            del_neighbors: IntMap::default(),
            three_way_shift_neighbors: cfg.three_way
                .then(ThreeWayNeighbors::default),
            four_way_shift_neighbors: cfg.four_way
                .then(FourWayNeighbors::default),
        };

        moves.init_del_neighbors();
        moves.init_loop_neighbors();
        moves
    }
}

impl<'a, T: LoopDecomposition, E: EnergyModel>
TryFrom<(&'a NucleotideVec, &T, &'a E, ShiftConfig)>
for AddDelShiftMoves<'a, E>
{
    type Error = String;

    fn try_from((sequence, pairings, model, cfg):
        (&'a NucleotideVec, &T, &'a E, ShiftConfig)
    ) -> Result<Self, Self::Error> {
        let ltab = LoopTable::try_from((sequence, pairings, model))?;
        Ok(AddDelShiftMoves::from((ltab, cfg)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Walker;
    use std::collections::HashSet;

    #[test]
    fn test_add_then_del_roundtrip() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GCUAACAACGGUCA");
        let pairings =       PairTable::try_from("..(.......)...").unwrap();

        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from((ltab, ShiftConfig {
            three_way: true,
            four_way: false,
        }));

        // Clone neighbor list so we don’t mutate while iterating
        let neighbors: Vec<(Move, i32)> = adm.propose_moves().collect();

        for (mv, de) in neighbors {
            let en0 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv, de, en0);

            // add pair
            let _ = adm.apply_move(&mv);
            let en1 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv.inverse(), -de, en1);

            let _ = adm.apply_move(&mv.inverse());
            let en2 = adm.current_energy();
            println!("({:?} {}) at energy: {}", mv.inverse().inverse(), de, en2);

            assert_eq!(en0, en2, "roundtrip energy mismatch");
        }
    }

    #[test]
    fn test_development_bug01() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GCAUAGCCCA");
        let pairings =       PairTable::try_from("..........").unwrap();

        let exp_nb1: Vec<(Move, i32)> = 
            vec![(Move::Add { i: 0, j: 6 }, 380),
                 (Move::Add { i: 0, j: 7 }, 390),
                 (Move::Add { i: 0, j: 8 }, 360),
                 (Move::Add { i: 1, j: 5 }, 430),
                 (Move::Add { i: 3, j: 9 }, 560)];

        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from((ltab, ShiftConfig {
            three_way: true,
            four_way: false,
        }));
        let en1 = adm.current_energy();
        let nb1: Vec<_> = adm.propose_moves().collect();
        assert_eq!(exp_nb1, nb1);

        let (del, add) = adm.apply_move(&Move::Add { i: 0, j: 7 });
        println!("Applied add");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::ShiftIK { i: 0, j: 7, k: 6 });
        println!("Applied shift");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::Del { i: 0, j: 6 });
        println!("Applied del");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }
        assert_eq!(en1, adm.current_energy());
        assert_eq!(exp_nb1, adm.propose_moves().collect::<Vec<_>>());
    }

    #[test]
    fn test_development_bug02() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GUACACGUCG");
        let pairings =       PairTable::try_from("..........").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from((ltab, ShiftConfig {
            three_way: true,
            four_way: false,
        }));
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

        let (del, add) = adm.apply_move(&Move::ShiftIK { i: 0, j: 8, k: 5 });
        println!("Applied shift");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::Del { i: 0, j: 5 });
        println!("Applied del");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }
        assert_eq!(en1, adm.current_energy());
    }


    #[test]
    fn test_development_bug03() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("GACGCUAUGU");
        let pairings =       PairTable::try_from("...(.....)").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from((ltab, ShiftConfig {
            three_way: true,
            four_way: false,
        }));
        let nb1: Vec<_> = adm.propose_moves().collect();

        for &(mv, d) in &nb1 {
            println!("pp {:?} {}", mv, d);
        }

        let pp1 = Move::Add { i: 4, j: 8 };
        let pp2 = Move::Del { i: 3, j: 9 };
        let pp3 = Move::ShiftIK { i: 3, j: 9, k: 7 };
        let pp4 = Move::ShiftJK { i: 3, j: 9, k: 0 };
        let pp5 = Move::ShiftJK { i: 3, j: 9, k: 1 };

        let ad1 = Move::Add { i: 1, j: 5 };
        let ad2 = Move::Add { i: 1, j: 7 };
        let ad3 = Move::Add { i: 2, j: 8 };
        let ad4 = Move::Add { i: 3, j: 7 };
        let ad5 = Move::ShiftJK { i: 0, j: 9, k: 1 };
        let ad6 = Move::ShiftJK { i: 0, j: 9, k: 3 };
        let ad7 = Move::ShiftIK { i: 0, j: 9, k: 4 };
        let ad8 = Move::ShiftIK { i: 0, j: 9, k: 5 };
        let ad9 = Move::ShiftIK { i: 0, j: 9, k: 7 };
        let ad0 = Move::Del { i: 0, j: 9 };

        let expected: HashSet<_> = [pp1, pp2, pp3, pp4, pp5].into_iter().collect();
        let actual: HashSet<_> = nb1.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);

        let (del, add) = adm.apply_move(&Move::ShiftJK{ i: 3, j: 9, k: 0 });

        println!("Applied shift");
        for &(mv, d) in &del {
            println!("rm {:?} {}", mv, d);
        }
        let expected: HashSet<_> = [pp2, pp3, pp4, pp5].into_iter().collect();
        let actual: HashSet<_> = del.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
        for &(mv, d) in &add {
            println!("up {:?} {}", mv, d);
        }
        let expected: HashSet<_> = [pp1, ad1, ad2, ad3, ad4, ad5, ad6, ad7, ad8, ad9, ad0].into_iter().collect();
        let actual: HashSet<_> = add.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_development_bug04() {
        let model = ViennaRNA::default();
        let sequence = NucleotideVec::from_lossy("AAAGCAAAAGCAAAAAAGAAAC");
        let pairings =       PairTable::try_from("...((....))......(...)").unwrap();
        let ltab = LoopTable::try_from((&sequence, &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from((ltab, ShiftConfig {
            three_way: true,
            four_way: true,
        }));
        let nb1: Vec<_> = adm.propose_moves().collect();
        let pp1 = Move::Del { i: 4, j: 9 };
        let pp2 = Move::Del { i: 17, j: 21 };
        let pp3 = Move::Del { i: 3, j: 10 };
        let pp4 = Move::ShiftILJK { i: 3, j: 10, k: 17, l: 21 };

        let am1 = Move::Del { i: 3, j: 21 };
        let am2 = Move::Del { i: 10, j: 17 };
        let am3 = Move::ShiftIKLJ { i: 3, j: 21, k: 10, l: 17 };

        for &(mv, d) in &nb1 {
            println!("pp {:?} {}", mv, d);
        }
        let expected: HashSet<_> = [pp1, pp2, pp3, pp4].into_iter().collect();
        let actual: HashSet<_> = nb1.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);

        let (del, add) = adm.apply_move(&pp4);
        println!("Applied shift");
        for &(mv, d) in &del {
            println!("rm {:?} {}", mv, d);
        }
        let expected: HashSet<_> = [pp2, pp3, pp4].into_iter().collect();
        let actual: HashSet<_> = del.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
        for &(mv, d) in &add {
            println!("up {:?} {}", mv, d);
        }
        let expected: HashSet<_> = [pp1, am1, am2, am3].into_iter().collect();
        let actual: HashSet<_> = add.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
    }


}

