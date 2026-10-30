use std::fmt;
use std::sync::Arc;
use nohash_hasher::IntMap;

use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_structure::MultiPairTable;

use crate::MoveSetError;
use crate::Move;
use crate::Neighbors;
use crate::movesets::loop_table::LoopTable;
use crate::movesets::three_way_shifts::ThreeWayNeighbors;
use crate::movesets::four_way_shifts::FourWayNeighbors;
use crate::ShiftMoves;


pub struct LoopNeighbors<E: EnergyModel> {
    loop_table: LoopTable<E>,
    add_neighbors: IntMap<usize, Neighbors>,
    del_neighbors: IntMap<usize, i32>, 
    three_way_shift_neighbors: ThreeWayNeighbors,
    four_way_shift_neighbors: FourWayNeighbors,
    shifts: ShiftMoves,
}

impl<E: EnergyModel> LoopNeighbors<E> {
    pub fn loop_table(&self) -> &LoopTable<E> {
        &self.loop_table
    }

    pub(crate) fn add_neighbors(&self) -> &IntMap<usize, Neighbors> {
        &self.add_neighbors
    }

    pub(crate) fn del_neighbors(&self) -> &IntMap<usize, i32> {
        &self.del_neighbors
    }

    pub(crate) fn three_way_shift_neighbors(&self) -> &ThreeWayNeighbors {
        &self.three_way_shift_neighbors
    }

    pub(crate) fn four_way_shift_neighbors(&self) -> &FourWayNeighbors {
        &self.four_way_shift_neighbors
    }

    /// Activation energy -> for add moves it is delta-E
    /// We should be able to change this for alternative 
    /// rate models, but do not break detailed balance.
    /// None if a new loop does not fit (see `LoopTable::fits`).
    fn get_add_activation_energy(&self, 
        index: usize,
        i: usize, 
        j: usize,
    ) -> Option<i32> {
        let ltab = &self.loop_table;
        let (combo, combo_energy) = ltab.get(index);
        let (outer, inner) = combo.split_loop(i, j);
        if !ltab.fits(&outer) || !ltab.fits(&inner) {
            return None;
        }
        let outer_energy = ltab.energy_of_loop(&outer);
        let inner_energy = ltab.energy_of_loop(&inner);
        Some((outer_energy + inner_energy) - combo_energy)
    }

    /// Returns how the free energy changes if the move is applied.
    fn get_del_activation_energy(&self, 
        i: usize, 
        j: usize,
    ) -> i32 {
        let (outer, outer_en) = self.loop_table.geti(i);
        let (inner, inner_en) = self.loop_table.geti(j);
        let combo = outer.join_loop(inner);
        let combo_energy = self.loop_table.energy_of_loop(&combo);
        combo_energy - (outer_en + inner_en)
    }

    fn init_del_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for (&i, &j) in ltab.pairs() {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
        }
    }

    #[inline(always)]
    fn maybe_remove_three_way(&mut self, idx: &usize) -> Neighbors {
        if self.shifts.three_way {
            self.three_way_shift_neighbors.remove(idx)
                .expect("Missing three-way shift neighbors.")
        } else {
            Vec::new()
        }
    }

    #[inline(always)]
    fn maybe_remove_four_way(&mut self, idx: &usize) -> Neighbors {
        if self.shifts.four_way {
            self.four_way_shift_neighbors.remove(idx)
                .expect("Missing four-way shift neighbors.")
        } else {
            Vec::new()
        }
    }

    #[inline(always)]
    fn maybe_three_way(&mut self, idx: usize) -> Neighbors {
        if self.shifts.three_way {
            self.three_way_shift_neighbors
                .compute_neighbors(&self.loop_table, idx)
                .clone()
        } else {
            Vec::new()
        }
    }

    #[inline(always)]
    fn maybe_four_way(&mut self, idx: usize) -> Neighbors {
        if self.shifts.four_way {
            self.four_way_shift_neighbors
                .compute_neighbors(&self.loop_table, idx)
                .clone()
        } else {
            Vec::new()
        }
    }

    fn init_loop_neighbors(&mut self) {
        for lli in 0..self.loop_table.loops_len() {
            self.get_add_neighbors_per_loop(lli);

            if self.shifts.three_way {
                let ltab = &self.loop_table;
                self.three_way_shift_neighbors
                    .compute_neighbors(ltab, lli);
            }

            if self.shifts.four_way {
                let ltab = &self.loop_table;
                self.four_way_shift_neighbors
                    .compute_neighbors(ltab, lli);
            }
        }
    }

    fn get_add_neighbors_per_loop(&mut self, index: usize) -> &Neighbors {
        let (combo, _) = self.loop_table.get(index);
        let unpaired: Vec<usize> = combo.stretches().flat_map(|s| s.unpaired).collect();
        let mut neighbors = Vec::new(); 
        for (idx_i, &i) in unpaired.iter().enumerate() {
            for &j in &unpaired[idx_i + 1..] {
                if j <= i + self.loop_table.min_hairpin_size() {
                    continue;
                }
                if self.loop_table.can_pair(i, j)
                    && let Some(barrier) = self.get_add_activation_energy(index, i, j)
                {
                    neighbors.push((Move::Add { i, j }, barrier)); 
                }
            }
        }
        self.add_neighbors.insert(index, neighbors);
        &self.add_neighbors[&index]
    }

    /// The delete moves of the pairs, each followed by the shift moves of
    /// the loop on its other side, if that loop is not among `changed`.
    fn pair_moves(&mut self, pairs: &[(usize, usize)], changed: &[usize], moves: &mut Neighbors) {
        for &(i, j) in pairs {
            let delta = self.get_del_activation_energy(i, j);
            self.del_neighbors.insert(i, delta);
            moves.push((Move::Del { i, j }, delta));
            if self.shifts.three_way || self.shifts.four_way {
                for side in [self.loop_table.loop_lookup(i), self.loop_table.loop_lookup(j)] {
                    if !changed.contains(&side) {
                        moves.extend(self.maybe_three_way(side));
                        moves.extend(self.maybe_four_way(side));
                    }
                }
            }
        }
    }

    /// The moves of new loops: their add moves, then their three-way and
    /// their four-way shifts.
    fn loop_moves(&mut self, loops: &[usize], moves: &mut Neighbors) {
        for &l in loops {
            moves.extend(self.get_add_neighbors_per_loop(l).iter().copied());
        }
        for &l in loops {
            moves.extend(self.maybe_three_way(l));
        }
        for &l in loops {
            moves.extend(self.maybe_four_way(l));
        }
    }

    pub(crate) fn apply_ext_move(&mut self) -> (Neighbors, Neighbors) {
        let index = self.loop_table.loop_lookup(0);
        let (old, _) = self.loop_table.get(index);
        let new = old.extended();
        let new_en = self.loop_table.energy_of_loop(&new);
        let new_pairs = new.index_pairs();
        let _ = self.loop_table.insert_loopentry(Some(index), (new, new_en));
        self.loop_table.extend_lookup(index);

        let mut new_moves = Vec::new();
        self.pair_moves(&new_pairs, &[index], &mut new_moves);
        // The extended loop: its energy changed, and the new base can pair.
        // No old move becomes invalid.
        self.loop_moves(&[index], &mut new_moves);
        (Vec::new(), new_moves)
    }

    pub(crate) fn apply_del_move(&mut self, i: usize, j: usize) -> (Neighbors, Neighbors) 
    {
        let o_index = self.loop_table.loop_lookup(i);
        let i_index = self.loop_table.loop_lookup(j);

        // Remove all old neighbors related to deletion.
        let delta = self.del_neighbors.remove(&i)
            .expect("Missing old del move.");
        let _ = self.add_neighbors.remove(&i_index)
            .expect("Missing old add moves from stale index");
        let old_tw_shift_outer = self.maybe_remove_three_way(&o_index);
        let old_tw_shift_inner = self.maybe_remove_three_way(&i_index);
        let old_fw_shift_outer = self.maybe_remove_four_way(&o_index);
        let old_fw_shift_inner = self.maybe_remove_four_way(&i_index);

        // Filter for all moves that cannot be regenerated. That is, 
        // they depend on the presence of the removed pair.
        let old_moves: Neighbors =
            [old_tw_shift_outer,
             old_tw_shift_inner,
             old_fw_shift_outer,
             old_fw_shift_inner].into_iter()
            .flatten()
            .filter(move |(mv, _)| mv.conflicts((i, j)))
            .chain(std::iter::once((Move::Del { i, j }, delta)))
            .collect();

        // Update the LoopTable.
        let (outer, o_en) = self.loop_table.get(o_index);
        let (inner, i_en) = self.loop_table.get(i_index);
        let combo = outer.join_loop(inner);
        let combo_pairs = combo.index_pairs();
        let c_en = o_en + i_en + delta;
        debug_assert_eq!(c_en, self.loop_table.energy_of_loop(&combo));

        let c_index = self.loop_table
            .insert_loopentry(Some(o_index), (combo, c_en));
        self.loop_table.mark_stale(i_index);
        self.loop_table.update_lookup(c_index);
        if j != self.loop_table.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }

        let mut new_moves = Vec::new();
        self.pair_moves(&combo_pairs, &[c_index], &mut new_moves);
        self.loop_moves(&[c_index], &mut new_moves);
        (old_moves, new_moves)
    }

    pub(crate) fn apply_add_move(&mut self, i: usize, j: usize) -> (Neighbors, Neighbors) 
    {
        // Get the original "combo" loop index
        let c_index = self.loop_table.loop_lookup(i);
        debug_assert_eq!(c_index, self.loop_table.loop_lookup(j), 
            "Missing loop_lookup entry for j.");

        // Remove all old neighbors related to addition.
        let old_add_moves = self.add_neighbors.remove(&c_index)
            .expect("Missing old add moves.");
        let old_tw_shift = self.maybe_remove_three_way(&c_index);
        let old_fw_shift = self.maybe_remove_four_way(&c_index);
        // Filter for all moves that cannot be regenerated. That is, 
        // they cross the newly added pair.
        let old_moves: Neighbors =
            [old_add_moves,
             old_tw_shift,
             old_fw_shift].into_iter()
            .flatten()
            .filter(move |(mv, _)| mv.conflicts((i, j)))
            .collect();

        // Update the LoopTable.
        let (combo, c_en) = self.loop_table.get(c_index);
        let combo_pairs = combo.index_pairs();
        let (outer, inner) = combo.split_loop(i, j);
        let o_en = self.loop_table.energy_of_loop(&outer);
        let i_en = self.loop_table.energy_of_loop(&inner);
        let delta = (o_en + i_en) - c_en; 

        self.loop_table.insert_pair(i, j);
        let o_index = self.loop_table.insert_loopentry(Some(c_index), (outer, o_en));
        let i_index = self.loop_table.insert_loopentry(None, (inner, i_en));
        self.loop_table.update_lookup(o_index);
        self.loop_table.update_lookup(i_index);

        let mut new_moves = Vec::new();
        self.pair_moves(&combo_pairs, &[o_index, i_index], &mut new_moves);
        self.del_neighbors.insert(i, -delta);
        new_moves.push((Move::Del { i, j }, -delta));
        self.loop_moves(&[o_index, i_index], &mut new_moves);
        (old_moves, new_moves)
    }

    pub(crate) fn apply_three_way_shift_move(
        &mut self, 
        mv: &Move, 
        k: usize,
    ) -> (Neighbors, Neighbors) {

        let (i, j) = mv.deleted_pair();
        let (p, q) = mv.added_pair();

        let o_index = self.loop_table.loop_lookup(i);
        let i_index = self.loop_table.loop_lookup(j);
        let k_index = self.loop_table.loop_lookup(k);
        debug_assert!(o_index == k_index || i_index == k_index);

        // Remove all old neighbors related to shift.
        let delta = self.del_neighbors.remove(&i)
            .expect("Missing del neighbor.");
        let old_add_init = self.add_neighbors.remove(&k_index)
            .expect("Missing old add moves from k_index.");
        let old_tw_shift_outer = self.three_way_shift_neighbors.remove(&o_index)
                .expect("Missing three-way shift neighbors.");
        let old_tw_shift_inner =  self.three_way_shift_neighbors.remove(&i_index)
                .expect("Missing three-way shift neighbors.");
        let old_fw_shift_outer = self.maybe_remove_four_way(&o_index);
        let old_fw_shift_inner = self.maybe_remove_four_way(&i_index);

        // Filter for all moves that cannot be regenerated. That is, 
        // they depend on the presence of the removed pair.
        let old_moves: Neighbors =
            old_add_init.into_iter()
            .filter(|(mv, _)| mv.conflicts((p, q)))
            .chain(old_tw_shift_outer) 
            .chain(old_tw_shift_inner)
            .chain(old_fw_shift_outer)
            .chain(old_fw_shift_inner)
            .filter(|(mv, _)| mv.conflicts((p, q)) || mv.conflicts((i, j)))
            .chain(std::iter::once((Move::Del { i, j }, delta)))
            .collect();

        let (outer, o_en) = self.loop_table.get(o_index);
        let (inner, i_en) = self.loop_table.get(i_index);

        let combo = outer.join_loop(inner);
        let combo_pairs = combo.index_pairs();
        let (new_outer, new_inner) = combo.split_loop(p, q);
        let new_o_en = self.loop_table.energy_of_loop(&new_outer);
        let new_i_en = self.loop_table.energy_of_loop(&new_inner);

        // This is a bit annoying for del_delta. 
        let c_en = o_en + i_en + delta;
        debug_assert_eq!(c_en, self.loop_table.energy_of_loop(&combo));
        let del_delta = c_en - (new_o_en + new_i_en);

        let _ = self.loop_table.insert_loopentry(Some(o_index), (new_outer, new_o_en));
        let _ = self.loop_table.insert_loopentry(Some(i_index), (new_inner, new_i_en));
        self.loop_table.update_lookup(o_index);
        self.loop_table.update_lookup(i_index);
        if j != self.loop_table.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        self.loop_table.insert_pair(p, q);

        let mut new_moves = Vec::new();
        self.pair_moves(&combo_pairs, &[o_index, i_index], &mut new_moves);
        self.del_neighbors.insert(p, del_delta);
        new_moves.push((Move::Del { i: p, j: q }, del_delta));
        self.loop_moves(&[o_index, i_index], &mut new_moves);
        (old_moves, new_moves)
    }

    pub(crate) fn apply_four_way_shift_move(
        &mut self, 
        mv: &Move, 
    ) -> (Neighbors, Neighbors) {
        let ((i, j), (k, l)) = mv.deleted_pairs();

        let (center_idx, merge1_idx, merge2_idx, inside, outer1, outer2) = 
            self.four_way_shift_neighbors.get_loops(&self.loop_table, mv);

        // Remove all old neighbors related to shift.
        let deld1 = self.del_neighbors.remove(&i).expect("Missing del neighbor.");
        let deld2 = self.del_neighbors.remove(&k).expect("Missing del neighbor.");
        let old_add_center = self.add_neighbors.remove(&center_idx)
            .expect("Missing old add moves");
        let old_add_merge1 = self.add_neighbors.remove(&merge1_idx)
            .expect("Missing old add moves");
        let old_add_merge2 = self.add_neighbors.remove(&merge2_idx)
            .expect("Missing old add moves");
        let old_tw_shift_center = self.maybe_remove_three_way(&center_idx);
        let old_tw_shift_merge1 = self.maybe_remove_three_way(&merge1_idx);
        let old_tw_shift_merge2 = self.maybe_remove_three_way(&merge2_idx);
        let old_fw_shift_center = self.maybe_remove_four_way(&center_idx); // definitely
        let old_fw_shift_merge1 = self.maybe_remove_four_way(&merge1_idx); // definitely
        let old_fw_shift_merge2 = self.maybe_remove_four_way(&merge2_idx); // definitely

        // Filter for all moves that cannot be regenerated. That is, 
        // they depend on the presence of the added pairs.
        let ((p, q), (m, n)) = mv.added_pairs();
        let old_moves: Neighbors = old_add_center.into_iter()
            .chain(old_add_merge1) 
            .chain(old_add_merge2) 
            .chain(old_tw_shift_center) 
            .chain(old_tw_shift_merge1) 
            .chain(old_tw_shift_merge2) 
            .chain(old_fw_shift_center) 
            .chain(old_fw_shift_merge1) 
            .chain(old_fw_shift_merge2) 
            .filter(|(mv, _)| mv.conflicts((p, q)) || mv.conflicts((m, n)))
            .chain(std::iter::once((Move::Del { i, j }, deld1)))
            .chain(std::iter::once((Move::Del { i: k, j: l }, deld2)))
            .collect();

        let inside_pairs = inside.index_pairs();
        let outer1_pairs = outer1.index_pairs();
        let outer2_pairs = outer2.index_pairs();
        let inside_en = self.loop_table.energy_of_loop(&inside);  
        let outer1_en = self.loop_table.energy_of_loop(&outer1);
        let outer2_en = self.loop_table.energy_of_loop(&outer2);

        let in_idx = self.loop_table.insert_loopentry(Some(center_idx), (inside, inside_en));
        let o1_idx = self.loop_table.insert_loopentry(Some(merge1_idx), (outer1, outer1_en));
        let o2_idx = self.loop_table.insert_loopentry(Some(merge2_idx), (outer2, outer2_en));
        self.loop_table.update_lookup(in_idx);
        self.loop_table.update_lookup(o1_idx);
        self.loop_table.update_lookup(o2_idx);
        if j != self.loop_table.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        if l != self.loop_table.delete_pair(&k) {
            panic!("Inconsistent pair-list entry.");
        }
        self.loop_table.insert_pair(p, q);
        self.loop_table.insert_pair(m, n);

        // The pairs of the new loops, the two new pairs once (in `inside`).
        let pairs: Vec<(usize, usize)> = inside_pairs.into_iter()
            .chain(outer1_pairs.into_iter().chain(outer2_pairs).filter(|&(a, _)| a != p && a != m))
            .collect();
        let mut new_moves = Vec::new();
        self.pair_moves(&pairs, &[in_idx, o1_idx, o2_idx], &mut new_moves);
        self.loop_moves(&[in_idx, o1_idx, o2_idx], &mut new_moves);
        (old_moves, new_moves)
    }
}

impl<E: EnergyModel> Clone for LoopNeighbors<E> {
    fn clone(&self) -> Self {
        Self {
            loop_table: self.loop_table.clone(),
            add_neighbors: self.add_neighbors.clone(),
            del_neighbors: self.del_neighbors.clone(),
            three_way_shift_neighbors: self.three_way_shift_neighbors.clone(),
            four_way_shift_neighbors: self.four_way_shift_neighbors.clone(),
            shifts: self.shifts,
        }
    }
}

impl<E: EnergyModel> fmt::Display for LoopNeighbors<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loop_table)
    }
}

impl<E: EnergyModel> From<(LoopTable<E>, ShiftMoves)> for LoopNeighbors<E> {
    fn from((loop_table, shifts): (LoopTable<E>, ShiftMoves)) -> Self {
        let mut moves = LoopNeighbors {
            loop_table,
            add_neighbors: IntMap::default(),
            del_neighbors: IntMap::default(),
            three_way_shift_neighbors: ThreeWayNeighbors::default(),
            four_way_shift_neighbors: FourWayNeighbors::default(),
            shifts,
        };

        moves.init_del_neighbors();
        moves.init_loop_neighbors();
        moves
    }
}

impl<E: EnergyModel> TryFrom<(NucleotideVec, &MultiPairTable, Arc<E>, ShiftMoves)>
for LoopNeighbors<E>
{
    type Error = MoveSetError;

    fn try_from((sequence, pairings, model, shifts):
        (NucleotideVec, &MultiPairTable, Arc<E>, ShiftMoves)
    ) -> Result<Self, Self::Error> {
        let ltab = LoopTable::try_from((sequence, pairings, model))?;
        Ok(LoopNeighbors::from((ltab, shifts)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::MultiPairTable;
    use ff_energy::Fuzzyfold;
    use ff_energy::NucleotideVec;
    use ff_energy::parameters::RNA_TURNER_2004;
    use std::collections::HashSet;
    use crate::MoveSet;

    macro_rules! setup_loop_table {
        ($name:ident, $seq:expr, $db:expr) => {
            let model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
            let sequence = NucleotideVec::try_from($seq).unwrap();
            let pairings = MultiPairTable::try_from($db)
                .expect("Invalid structure");
            let $name = LoopTable::try_from((sequence, &pairings, Arc::new(model)))
                .expect("Invalid sequence/structure combination");
        };
    }

    #[test]
    fn test_add_then_del_roundtrip() {
        setup_loop_table!(ltab, "GCUAACAACGGUCA", 
                                "..(.......)...");
        let mut adm = LoopNeighbors::from((ltab, ShiftMoves { three_way: true, four_way: false }));

        // Collect the moves first: applying them changes the move set.
        let neighbors: Vec<(Move, i32)> = adm.propose_moves().collect();

        for (mv, de) in neighbors {
            let en0 = adm.current_energy();
            let _ = adm.apply_move(&mv);
            // The delta of a shift move is its barrier, not the energy change.
            if matches!(mv, Move::Add { .. } | Move::Del { .. }) {
                assert_eq!(adm.current_energy(), en0 + de, "{mv:?}");
            }
            let _ = adm.apply_move(&mv.inverse());
            assert_eq!(adm.current_energy(), en0, "roundtrip energy mismatch");
        }
    }

    #[test]
    fn test_development_bug01() {
        setup_loop_table!(ltab, "GCAUAGCCCA", 
                                "..........");
        let exp_nb1: Vec<(Move, i32)> = 
            vec![(Move::Add { i: 0, j: 6 }, 380),
                 (Move::Add { i: 0, j: 7 }, 390),
                 (Move::Add { i: 0, j: 8 }, 360),
                 (Move::Add { i: 1, j: 5 }, 430),
                 (Move::Add { i: 3, j: 9 }, 560)];

        let mut adm = LoopNeighbors::from((ltab, ShiftMoves { three_way: true, four_way: false }));
        let en1 = adm.current_energy();
        let nb1: Vec<_> = adm.propose_moves().collect();
        assert_eq!(exp_nb1, nb1);

        let _ = adm.apply_move(&Move::Add { i: 0, j: 7 });
        let _ = adm.apply_move(&Move::ShiftIK { i: 0, j: 7, k: 6 });
        let _ = adm.apply_move(&Move::Del { i: 0, j: 6 });
        assert_eq!(en1, adm.current_energy());
        assert_eq!(exp_nb1, adm.propose_moves().collect::<Vec<_>>());
    }

    #[test]
    fn test_development_bug02() {
        setup_loop_table!(ltab, "GUACACGUCG", 
                                "..........");
        let mut adm = LoopNeighbors::from((ltab, ShiftMoves { three_way: true, four_way: false }));
        let en1 = adm.current_energy();
        let nb1: Vec<_> = adm.propose_moves().collect();

        let _ = adm.apply_move(&Move::Add { i: 0, j: 8 });
        let _ = adm.apply_move(&Move::ShiftIK { i: 0, j: 8, k: 5 });
        let _ = adm.apply_move(&Move::Del { i: 0, j: 5 });
        assert_eq!(en1, adm.current_energy());
        let mut after: Vec<_> = adm.propose_moves().collect();
        let mut before = nb1;
        after.sort_by_key(|&(mv, _)| format!("{mv:?}"));
        before.sort_by_key(|&(mv, _)| format!("{mv:?}"));
        assert_eq!(before, after);
    }


    #[test]
    fn test_development_bug03() {
        let model = Fuzzyfold::default();
        let sequence = NucleotideVec::try_from("GACGCUAUGU").unwrap();
        let pairings =       MultiPairTable::try_from("...(.....)").unwrap();
        let ltab = LoopTable::try_from((sequence, &pairings, Arc::new(model))).unwrap();
        let mut adm = LoopNeighbors::from((ltab, ShiftMoves { three_way: true, four_way: false }));
        let nb1: Vec<_> = adm.propose_moves().collect();

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

        let expected: HashSet<_> = [pp2, pp3, pp4, pp5].into_iter().collect();
        let actual: HashSet<_> = del.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
        let expected: HashSet<_> = [pp1, ad1, ad2, ad3, ad4, ad5, ad6, ad7, ad8, ad9, ad0].into_iter().collect();
        let actual: HashSet<_> = add.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_development_bug04() {
        let model = Fuzzyfold::default();
        let sequence = NucleotideVec::try_from_rna("AAAGCAAAAGCAAAAAAGAAAC").unwrap();
        let pairings =       MultiPairTable::try_from("...((....))......(...)").unwrap();
        let ltab = LoopTable::try_from((sequence, &pairings, Arc::new(model))).unwrap();
        let mut adm = LoopNeighbors::from((ltab, ShiftMoves::ALL));
        let nb1: Vec<_> = adm.propose_moves().collect();
        let pp1 = Move::Del { i: 4, j: 9 };
        let pp2 = Move::Del { i: 17, j: 21 };
        let pp3 = Move::Del { i: 3, j: 10 };
        let pp4 = Move::ShiftILJK { i: 3, j: 10, k: 17, l: 21 };

        let am1 = Move::Del { i: 3, j: 21 };
        let am2 = Move::Del { i: 10, j: 17 };
        let am3 = Move::ShiftIKLJ { i: 3, j: 21, k: 10, l: 17 };

        let expected: HashSet<_> = [pp1, pp2, pp3, pp4].into_iter().collect();
        let actual: HashSet<_> = nb1.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);

        let (del, add) = adm.apply_move(&pp4);
        let expected: HashSet<_> = [pp2, pp3, pp4].into_iter().collect();
        let actual: HashSet<_> = del.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
        let expected: HashSet<_> = [pp1, am1, am2, am3].into_iter().collect();
        let actual: HashSet<_> = add.iter().map(|(mv, _)| mv).cloned().collect();
        assert_eq!(actual, expected);
    }


}

