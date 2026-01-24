use std::fmt;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_energy::LoopDecomposition;
use ff_energy::NearestNeighborLoop;
use crate::Move;
use crate::movesets::LoopTable;

type Pair = (NAIDX, NAIDX);
type Moves = Vec<(Move, i32)>;
type LoopEntry = (NearestNeighborLoop, i32);

pub struct AddDelShiftMoves<'a, E: EnergyModel> {
    loop_table: LoopTable<'a, E>,
    /// registry index to list of (i, j, deltaE)
    add_neighbors: IntMap<usize, Moves>,
    /// pair id to deltaE
    del_neighbors: IntMap<NAIDX, i32>, 
    /// registry index to list of (i, j, deltaE)
    shift_neighbors: IntMap<usize, Moves>,
}

impl<'a, E: EnergyModel> fmt::Display for AddDelShiftMoves<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.loop_table)
    }
}

impl<'a, E: EnergyModel> AddDelShiftMoves<'a, E> {
    pub fn loop_table(&self) -> &LoopTable<'a, E> {
        &self.loop_table
    }

    pub fn add_neighbors(&self) -> &IntMap<usize, Moves> {
        &self.add_neighbors
    }

    pub fn shift_neighbors(&self) -> &IntMap<usize, Moves> {
        &self.shift_neighbors
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
    ) -> Moves {
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
        for lli in 0..ltab.loops.len() {
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
 
    fn init_shift_neighbors(&mut self) {
        let ltab = &self.loop_table;
        for lli in 0..ltab.loops.len() {
            let neighbors = self.get_shift_neighbors_per_loop(lli);
            self.shift_neighbors.insert(lli, neighbors);
        }
    }

    fn get_shift_neighbors_per_loop(&self, index: usize) -> Moves {
        let ltab = &self.loop_table;
        let (combo, combo_energy) = ltab.get(index);

        let mut loopdict: IntMap<NAIDX, (Pair, LoopEntry)> = IntMap::default();
        let mut neighbors = Vec::new(); 

        match combo {
            //TODO early exit!
            NearestNeighborLoop::Hairpin { closing: (i, j) } => {
                self.shift_loops_insert(*i, *j, *i as usize, (combo, combo_energy), &mut loopdict);
                self.shift_iter(*i as usize + 1, *j as usize, &loopdict, &mut neighbors);
            }
            NearestNeighborLoop::Interior { closing: (i, j), inner: (p, q) } => {
                self.shift_loops_insert(*i, *j, *i as usize, (combo, combo_energy), &mut loopdict);
                self.shift_loops_insert(*p, *q, *q as usize, (combo, combo_energy), &mut loopdict);

                self.shift_iter(*i as usize + 1, *p as usize, &loopdict, &mut neighbors);
                self.shift_iter(*q as usize + 1, *j as usize, &loopdict, &mut neighbors);
            },
            NearestNeighborLoop::Multibranch { closing: (i, j), branches } => {
                self.shift_loops_insert(*i, *j, *i as usize, (combo, combo_energy), &mut loopdict);
                for &(p, q) in branches {
                    self.shift_loops_insert(p, q, q as usize, (combo, combo_energy), &mut loopdict);
                }

                let mut start = *i as usize;
                for &(p, q) in branches {
                    self.shift_iter(start + 1, p as usize, &loopdict, &mut neighbors);
                    start = q as usize;
                }
                self.shift_iter(start + 1, *j as usize, &loopdict, &mut neighbors);
            },
            NearestNeighborLoop::Exterior { ends: (p5, p3), branches } => {
                for &(p, q) in branches {
                    self.shift_loops_insert(p, q, q as usize, (combo, combo_energy), &mut loopdict);
                }
                let mut start = *p5 as usize;
                for &(p, q) in branches {
                    self.shift_iter(start, p as usize, &loopdict, &mut neighbors);
                    start = q as usize + 1;
                }
                if start != *p5 as usize {
                    self.shift_iter(start, *p3 as usize + 1, &loopdict, &mut neighbors);
                }
            },
            _ => {panic!("no shift move for loop type!")},
        }
        neighbors
    }

    fn shift_loops_insert(
        &self,
        i: NAIDX,
        j: NAIDX,
        l: usize,
        (center, center_energy): (&NearestNeighborLoop, &i32),
        loopdict: &mut IntMap<NAIDX, (Pair, LoopEntry)>,
    ) {
        let ltab = &self.loop_table;
        let (shift, shift_energy) = ltab.geti(l); 
        let combo = if l == i as usize {
            shift.join_loop(center) 
        } else { 
            center.join_loop(shift)
        };
        //let combo = shift.join_loop(center); 
        let combo_energy = center_energy + shift_energy;
        loopdict.insert(i, ((i, j), (combo.clone(), combo_energy)));
        loopdict.insert(j, ((i, j), (combo, combo_energy)));
    }

    fn shift_iter(
        &self,
        p5a: usize,
        p3: usize,
        loopdict: &IntMap<NAIDX, (Pair, LoopEntry)>,
        neighbors: &mut Moves,
    ) {
        let ltab = &self.loop_table;
        let p5 = p5a.checked_sub(1).unwrap_or(0);
        let u5 = p5 as NAIDX;
        let u3 = p3 as NAIDX;
        for k in p5a..p3 {
            for (&p, ((pi, pj), (shift_combo_loop, shift_combo_energy))) in loopdict.iter() {
                if (p == u5 && k <= p5 + ltab.min_hairpin_size()) || 
                    (p == u3 && k + ltab.min_hairpin_size() >= p3) {
                        continue;
                } else if ltab.can_pair(k, p as usize) {
                    let nk = k as NAIDX;
                    let (i, j) = if p < nk { (p, nk) } else { (nk, p) };
                    let delta = self.shift_delta(
                        shift_combo_loop, *shift_combo_energy, i, j
                    );
                    let mv = if p == *pi {
                        Move::ShiftIK { i: *pi, j: *pj, k: nk }
                    } else {
                        Move::ShiftJK { i: *pi, j: *pj, k: nk }
                    };
                    neighbors.push((mv, delta));
                }
            }
        }
    }

    #[inline(always)]
    fn shift_delta(&self,
        combo: &NearestNeighborLoop,
        combo_energy: i32,
        split_a: NAIDX,
        split_b: NAIDX,
    ) -> i32 {
        let (s_outer, s_inner) = combo.split_loop(split_a, split_b);
        let s_outer_energy = self.loop_table.energy_of_loop(&s_outer);
        let s_inner_energy = self.loop_table.energy_of_loop(&s_inner);
        (s_outer_energy + s_inner_energy) - combo_energy
    }

    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let delta = self.del_neighbors.remove(&i).expect("Missing pair_list entry.");

        let ltab = &mut self.loop_table;
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        let o_index = ltab.loop_lookup(i as usize);
        let i_index = ltab.loop_lookup(j as usize);

        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);
        let combo = outer.join_loop(inner);
        let combo_energy = (o_en + i_en) + delta;
        debug_assert_eq!(combo_energy, ltab.energy_of_loop(&combo));

        let cpairs = &combo.pairs();

        ltab.loops[o_index] = (combo, combo_energy);
        ltab.update_lookup(o_index);
        ltab.stale.push(i_index);
        ltab.energy += delta;

        let _ = self.add_neighbors.remove(&i_index).expect("at least empty list.");
        let old_outer_shift_moves = self.shift_neighbors.remove(&o_index).expect("Old combo moves");
        let old_inner_shift_moves = self.shift_neighbors.remove(&i_index).expect("Old combo moves");
        let mut old_moves: Moves = old_outer_shift_moves
            .into_iter()
            .chain(old_inner_shift_moves)
            .filter(|(smv, _)| smv.conflicts((i, j)))
            .collect();
        old_moves.push((Move::Del { i, j }, delta));


        // Those include the neighbors
        let update_add_moves = self.get_add_neighbors_per_loop(o_index);
        let update_shift_moves = self.get_shift_neighbors_per_loop(o_index);
        let cap = update_add_moves.len() + update_shift_moves.len();
        let mut new_moves = self.update_del_moves(cpairs, cap);
        new_moves.extend(
            update_add_moves.iter().cloned()
            .chain(update_shift_moves.iter().cloned()));

        self.add_neighbors.insert(o_index, update_add_moves);
        self.shift_neighbors.insert(o_index, update_shift_moves);

        (old_moves, new_moves)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX) -> (Moves, Moves) 
    {
        let ltab = &mut self.loop_table;

        let c_index = ltab.loop_lookup(i as usize);
        debug_assert_eq!(c_index, ltab.loop_lookup(j as usize), "Missing loop_lookup entry for j.");

        let (combo, c_en) = ltab.get(c_index).clone();
        let mut old_moves: Moves = self.add_neighbors.remove(&c_index).expect("Old combo moves")
            .into_iter()
            .filter(|(amv, _)| amv.conflicts((i, j)))
            .collect();
        let old_shift_moves = self.shift_neighbors.remove(&c_index).expect("Old shift moves");
        old_moves.extend(
            old_shift_moves
            .into_iter()
            .filter(|(smv, _)| smv.conflicts((i, j)))
            .collect::<Moves>()
        );

        let (outer, inner) = combo.split_loop(i, j);
        
        let o_en = ltab.energy_of_loop(&outer);
        ltab.loops[c_index] = (outer, o_en);
        let i_en = ltab.energy_of_loop(&inner);
        let i_index = ltab.insert_loop(inner, i_en);

        let delta = (o_en + i_en) - c_en;
        ltab.energy += delta;
        ltab.insert_pair(i, j);
        ltab.update_lookup(c_index);
        ltab.update_lookup(i_index);
        let outer_add_moves = self.get_add_neighbors_per_loop(c_index);
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index);
        let outer_shift_moves = self.get_shift_neighbors_per_loop(c_index);
        let inner_shift_moves = self.get_shift_neighbors_per_loop(i_index);

        let cap = outer_add_moves.len() 
            + inner_add_moves.len() 
            + outer_add_moves.len()
            + inner_shift_moves.len()
            + outer_shift_moves.len() + 1;
        let mut new_moves = self.update_del_moves(&combo.pairs(), cap);
        new_moves.push((Move::Del{ i, j }, -delta));
        new_moves.extend(
            outer_add_moves.iter().cloned()
            .chain(inner_add_moves.iter().cloned())
            .chain(outer_shift_moves.iter().cloned())
            .chain(inner_shift_moves.iter().cloned()));

        self.add_neighbors.insert(c_index, outer_add_moves);
        self.add_neighbors.insert(i_index, inner_add_moves);
        self.shift_neighbors.insert(c_index, outer_shift_moves);
        self.shift_neighbors.insert(i_index, inner_shift_moves);
        self.del_neighbors.insert(i, -delta);

        (old_moves, new_moves)
    }

    pub fn apply_shift_move(
        &mut self, 
        mv: &Move, 
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
        (p, q): (NAIDX, NAIDX)
    ) -> (Moves, Moves) {
        let ltab = &mut self.loop_table;

        let delta = self.del_neighbors.remove(&i).expect("Missing pair_list entry.");
        if j != ltab.delete_pair(&i) {
            panic!("Inconsistent pair-list entry.");
        }
        let o_index = ltab.loop_lookup(i as usize);
        let i_index = ltab.loop_lookup(j as usize);
        let k_index = ltab.loop_lookup(k as usize);
        debug_assert!(o_index == k_index || i_index == k_index);

        let mut old_moves: Vec<_> = self.add_neighbors.remove(&k_index).expect("Old kloop moves")
            .into_iter()
            .filter(|(amv, _)| amv.conflicts(mv.added_pair()))
            .collect();
        let old_inner_shift_moves = self.shift_neighbors.remove(&i_index).expect("Old shift moves");
        let old_outer_shift_moves = self.shift_neighbors.remove(&o_index).expect("Old shift moves");
        old_moves.extend(
            old_inner_shift_moves.into_iter()
            .chain(old_outer_shift_moves)
            .filter(|(smv, _)| smv.conflicts(mv.added_pair()) 
                || smv.conflicts(mv.deleted_pair()))
            .collect::<Vec<_>>());
        old_moves.push((Move::Del { i, j }, delta));

        let (outer, o_en) = ltab.get(o_index);
        let (inner, i_en) = ltab.get(i_index);
        let combo = outer.join_loop(inner);
        let combo_energy = (o_en + i_en) + delta;
        debug_assert_eq!(combo_energy, ltab.energy_of_loop(&combo));
        let cpairs = &combo.pairs();

        let (new_outer, new_inner) = combo.split_loop(p, q);
        let new_o_en = ltab.energy_of_loop(&new_outer);
        let new_i_en = ltab.energy_of_loop(&new_inner);
        ltab.energy += (new_o_en + new_i_en) - (o_en + i_en);
        ltab.loops[o_index] = (new_outer, new_o_en);
        ltab.loops[i_index] = (new_inner, new_i_en);
        ltab.insert_pair(p, q);
        ltab.update_lookup(o_index);
        ltab.update_lookup(i_index);


        let outer_add_moves = self.get_add_neighbors_per_loop(o_index);
        let inner_add_moves = self.get_add_neighbors_per_loop(i_index);
        let outer_shift_moves = self.get_shift_neighbors_per_loop(o_index);
        let inner_shift_moves = self.get_shift_neighbors_per_loop(i_index);
        let cap = outer_add_moves.len() 
            + inner_add_moves.len() 
            + outer_add_moves.len()
            + inner_shift_moves.len()
            + outer_shift_moves.len();
        let mut new_moves = self.update_del_moves(cpairs, cap);

        let deld = combo_energy - (new_o_en + new_i_en);
        new_moves.extend(
            outer_add_moves.iter().cloned()
            .chain(inner_add_moves.iter().cloned())
            .chain(outer_shift_moves.iter().cloned())
            .chain(inner_shift_moves.iter().cloned()));
        new_moves.push((Move::Del { i: p, j: q }, deld));


        self.add_neighbors.insert(o_index, outer_add_moves);
        self.add_neighbors.insert(i_index, inner_add_moves);
        self.shift_neighbors.insert(o_index, outer_shift_moves);
        self.shift_neighbors.insert(i_index, inner_shift_moves);
        self.del_neighbors.insert(p, deld);

        (old_moves, new_moves)
    }

}

impl<'a, E: EnergyModel> From<LoopTable<'a, E>> for AddDelShiftMoves<'a, E> {
    fn from(loop_table: LoopTable<'a, E>) -> Self {
        let mut moves = AddDelShiftMoves {
            loop_table,
            add_neighbors: IntMap::default(),
            del_neighbors: IntMap::default(),
            shift_neighbors: IntMap::default(),
        };
        moves.init_add_neighbors();
        moves.init_del_neighbors();
        moves.init_shift_neighbors();
        moves
    }
}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a [Base], &T, &'a E)> 
for AddDelShiftMoves<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a E)
    ) -> Result<Self, Self::Error> {
        Ok(AddDelShiftMoves::from(LoopTable::try_from((sequence, pairings, model))?))
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
        let sequence = NucleotideVec::from_lossy("GCUAACAACGGUCA");
        let pairings =       PairTable::try_from("..(.......)...").unwrap();

        let ltab = LoopTable::try_from((&sequence[..], &pairings, &model)).unwrap();
        let mut adm = AddDelShiftMoves::from(ltab);

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

        let mut adm = AddDelShiftMoves::try_from((&sequence[..], &pairings, &model)).unwrap();
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
        let mut adm = AddDelShiftMoves::try_from((&sequence[..], &pairings, &model)).unwrap();
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
        let mut adm = AddDelShiftMoves::try_from((&sequence[..], &pairings, &model)).unwrap();
        let nb1: Vec<_> = adm.propose_moves().collect();
        for (mv, d) in nb1 {
            println!("pp {:?} {}", mv, d);
        }

        let (del, add) = adm.apply_move(&Move::ShiftJK{ i: 3, j: 9, k: 0 });
        println!("Applied shift");
        for (mv, d) in del {
            println!("rm {:?} {}", mv, d);
        }
        for (mv, d) in add {
            println!("up {:?} {}", mv, d);
        }
    }


}

