use std::fmt;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::NearestNeighborLoop;
use ff_energy::LoopDecomposition;
use ff_energy::EnergyModel;
use ff_energy::Base;

type LoopEntry = (NearestNeighborLoop, i32);
type MoveEnergies = Vec<(NAIDX, NAIDX, i32)>;

struct LoopCache<'a, E: EnergyModel> {
    sequence: &'a [Base],
    model: &'a E,
    loops: Vec<LoopEntry>,
    stale: Vec<usize>,
}

impl<'a, E: EnergyModel> LoopCache<'a, E> {

    pub fn new(sequence: &'a [Base], model: &'a E) -> Self {
        Self { 
            sequence,
            model,
            loops: Vec::new(),
            stale: Vec::new(),
        }
    }

    #[inline]
    pub fn get(&self, idx: usize) -> &LoopEntry {
        self.loops
            .get(idx)
            .expect("Invalid LoopCache index")
    }

    pub fn update(&mut self, idx: usize, new: NearestNeighborLoop) -> (usize, i32) {
        let energy = self.model.energy_of_loop(self.sequence, &new);
        self.loops[idx] = Some((new, energy));
        (idx, energy)
    }

    pub fn insert_loop(&mut self, combo: NearestNeighborLoop) -> (usize, i32) {
        let energy = self.model.energy_of_loop(self.sequence, &combo);

        if let Some(idx) = self.stale.pop() {
            self.loops[idx] = (combo, energy);
            (idx, energy)
        } else {
            let idx = self.loops.len();
            self.loops.push((combo, energy));
            (idx, energy)
        }
    }

    pub fn calc_pair_energy(&self, outer_index: usize, inner_index: usize) -> i32 {
        let (outer, o_en) = self.get(outer_index);
        let (inner, i_en) = self.get(inner_index);
        let combo = outer.join_loop(inner);
        let combo_energy = self.model.energy_of_loop(self.sequence, &combo);
        combo_energy - (o_en + i_en)
    }

    pub fn apply_delete_move(&mut self, outer_index: usize, inner_index: usize, delta: i32) -> usize {
        let (outer, o_en) = self.get(outer_index);
        let (inner, i_en) = self.get(inner_index);
        let combo = outer.join_loop(inner);
        let combo_energy = (o_en + i_en) - delta;
        // re-use outer_index for the new loop.
        self.loops[outer_index] = (combo, combo_energy);
        self.stale.push(inner_index);
        outer_index
    }

    pub fn apply_addition_move(&mut self, 
        c_index: usize, 
        combo: NearestNeighborLoop,
        c_en: i32, 
        i: NAIDX, j: NAIDX
    ) -> (usize, usize, i32) {
        let (outer, inner) = combo.split_loop(i, j);

        //NOTE: could look delta up directly by searching loop_list.
        let o_en = self.model.energy_of_loop(self.sequence, &outer);
        self.loops[c_index] = (outer, o_en);
        let (i_index, i_en) = self.insert_loop(inner);

        // How does the energy change if we apply the base-pair move.
        let delta = (o_en + i_en) - c_en;
        (c_index, i_index, -delta)
    }

    fn get_loop_neighbors(&self, index: usize) -> MoveEnergies {
        let (combo, energy) = self.get(index);
        let unpaired = combo.unpaired_indices();

        let mut neighbors = Vec::new(); 
        for (idx_i, &i) in unpaired.iter().enumerate() {
            for &j in &unpaired[idx_i + 1..] {
                if j <= i + self.model.min_hairpin_size() {
                    continue;
                }
                if self.model.can_pair(self.sequence[i], self.sequence[j]) {
                    let (outer, inner) = combo.split_loop(i as NAIDX, j as NAIDX);
                    let outer_energy = self.model.energy_of_loop(self.sequence, &outer);
                    let inner_energy = self.model.energy_of_loop(self.sequence, &inner);
                    // How does the free energy change if the move is applied.
                    let delta = (outer_energy + inner_energy) - energy;
                    neighbors.push((i as NAIDX, j as NAIDX, delta));
                }
            }
        }
        neighbors
    }
}


pub struct LoopStructure<'a, E: EnergyModel> {
    registry: LoopCache<'a, E>,
    /// From sequence index to registry index.
    loop_lookup: Vec<Option<usize>>, 
    /// registry index to list of (i, j, deltaE)
    loop_neighbors: IntMap<usize, MoveEnergies>,
    /// Current pairs, i<j where i is the id.
    pair_list: IntMap<NAIDX, NAIDX>,
    /// pair id to deltaE
    pair_neighbors: IntMap<NAIDX, i32>, 
    /// free energy of structure
    energy: i32,
}

impl<'a, E: EnergyModel> Clone for LoopStructure<'a, E> {
    fn clone(&self) -> Self {
        Self {
            registry: LoopCache {
                sequence: self.registry.sequence,
                model: self.registry.model,
                loops: self.registry.loops.clone(),
                stale: self.registry.stale.clone(),
            },
            loop_lookup: self.loop_lookup.clone(),
            loop_neighbors: self.loop_neighbors.clone(),
            pair_list: self.pair_list.clone(),
            pair_neighbors: self.pair_neighbors.clone(),
            energy: self.energy,
        }
    }
}

impl<'a, E: EnergyModel> fmt::Debug for LoopStructure<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructure")
            .field("loop_lookup", &self.loop_lookup)
            .field("num_pairs", &self.pair_list.len())
            .field("num_add_neighbors", &self.loop_neighbors.values().map(|v| v.len()).sum::<usize>())
            .field("num_del_neighbors", &self.pair_neighbors.len())
            .finish()
    }
}

impl<'a, E: EnergyModel> fmt::Display for LoopStructure<'a, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Convert sequence to string
        let mut dbr = vec!['.'; self.registry.sequence.len()];
        for (i, j) in &self.pair_list {
            dbr[*i as usize] = '(';
            dbr[*j as usize] = ')';
        }
        let dbr_str: String = dbr.into_iter().collect();
        write!(f, "{}", dbr_str)
    }
}

impl<'a, E: EnergyModel> LoopStructure<'a, E> {
    /// Return all add neighbors, including an index that 
    /// is necessary to access the actual loop via loop_lookup.
    pub fn get_add_neighbors_per_loop(&self) -> &IntMap<usize, MoveEnergies> {
        &self.loop_neighbors
    }

    pub fn len(&self) -> usize {
        self.pair_neighbors.len() + self.loop_neighbors.len()
    }
 
    pub fn is_empty(&self) -> bool {
        self.pair_neighbors.is_empty() && self.loop_neighbors.is_empty()
    }

    /// Return all remove neighbors, where all i, j are also
    /// the indices to access the outer/inner loop via loop_lookup.
    pub fn get_del_neighbors(&self) -> MoveEnergies {
        self.pair_neighbors
            .iter()
            .map(|(&i, &delta_e)| (i, self.pair_list[&i], delta_e))
            .collect()
    }

    /// A pair-table like structure, where each position points to 
    /// exactly one loop. 
    pub fn loop_lookup(&self) -> &Vec<Option<usize>> {
        &self.loop_lookup
    }

    pub fn energy(&self) -> i32 {
        self.energy
    }

    fn update_pair_neighbors(&mut self,
        pairs: &[(NAIDX, NAIDX)]
    ) -> MoveEnergies
    {
        let mut change = Vec::new();
        for &(i, j) in pairs {
            let o_index = self.loop_lookup[i as usize].unwrap();
            let i_index = self.loop_lookup[j as usize].unwrap();
            let delta = self.registry.calc_pair_energy(o_index, i_index);
            self.pair_neighbors.insert(i, delta);
            change.push((i, j, delta));
        }
        change
    }
  
    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX
    ) -> (MoveEnergies, MoveEnergies) 
    {
        if j != self.pair_list.remove(&i).expect("Missing pair_list entry.") {
            panic!("Inconsistent pair-list entry.");
        }

        let o_index = self.loop_lookup[i as usize].unwrap();
        let i_index = self.loop_lookup[j as usize].unwrap();
        let delta = self.pair_neighbors.remove(&i).expect("Missing pair_neighbors entry."); 
        self.energy += delta;

        let c_index = self.registry.apply_delete_move(o_index, i_index, -delta);
        assert_eq!(c_index, o_index); // by convention.

        let loop_neighbors = self.registry.get_loop_neighbors(o_index);
        self.loop_neighbors.insert(o_index, loop_neighbors.clone());
        self.loop_neighbors.remove(&i_index).expect("at least empty list.");

        let (combo, _) = self.registry.get(o_index);
        for k in &combo.inclusive_unpaired_indices() {
            debug_assert!(self.loop_lookup[*k] == Some(o_index) || self.loop_lookup[*k] == Some(i_index));
            self.loop_lookup[*k] = Some(o_index);
        }

        let pair_changes = self.update_pair_neighbors(&combo.pairs());
        (loop_neighbors, pair_changes)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX
    ) -> (MoveEnergies, MoveEnergies, MoveEnergies) 
    {
        let c_index = self.loop_lookup[i as usize].unwrap();
        debug_assert_eq!(Some(c_index), self.loop_lookup[j as usize], "Missing loop_lookup entry for j.");
        let (combo, c_en) = self.registry.get(c_index);
        let combo_pairs = &combo.pairs();
        // How does the energy change if we apply the base-pair move.
        let (o_id, i_id, delta) = self.registry.apply_addition_move(
            c_index, combo.to_owned(), *c_en, i, j
        );
        self.energy -= delta;

        let new_outer_add_neighbors = self.registry.get_loop_neighbors(o_id);
        let new_inner_add_neighbors = self.registry.get_loop_neighbors(i_id);
        self.loop_neighbors.insert(o_id, new_outer_add_neighbors.clone());
        self.loop_neighbors.insert(i_id, new_inner_add_neighbors.clone());
        self.pair_list.insert(i, j);
        self.pair_neighbors.insert(i, delta);

        let (outer, _) = self.registry.get(o_id);
        for k in &outer.inclusive_unpaired_indices() {
            self.loop_lookup[*k] = Some(o_id);
        }
        let (inner, _) = self.registry.get(i_id);
        for k in &inner.inclusive_unpaired_indices() {
            self.loop_lookup[*k] = Some(i_id);
        }

        let mut pair_changes = self.update_pair_neighbors(combo_pairs);
        pair_changes.push((i, j, delta));

        (new_outer_add_neighbors, new_inner_add_neighbors, pair_changes)
    }

    pub fn apply_ext_move(&mut self) -> (MoveEnergies, MoveEnergies) {
        let index = self.loop_lookup[0usize].unwrap();
        let &(ref old, old_en) = self.registry.get(index);
        let new = match old {
            NearestNeighborLoop::Exterior { ends: (p5, p3), branches } => {
                let p3u = *p3 as usize;
                debug_assert!(self.loop_lookup[p3u + 1].is_none());
                self.loop_lookup[p3u + 1] = Some(index);
                NearestNeighborLoop::Exterior { ends: (*p5, p3 + 1), branches: branches.clone() }
            },
            _ => panic!("should have been exterior loop"),
        };
        let (_, new_en) = self.registry.update(index, new);
        self.energy += new_en - old_en;

        let loop_neighbors = self.registry.get_loop_neighbors(index);
        self.loop_neighbors.insert(index, loop_neighbors.clone());

        let (new, _) = self.registry.get(index);
        let pair_changes = self.update_pair_neighbors(&new.pairs());
        (loop_neighbors, pair_changes)
    }

}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a [Base], &T, &'a E)> for LoopStructure<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a E)
    ) -> Result<Self, Self::Error> {
        let mut registry = LoopCache::new(sequence, model);
        let mut loop_lookup: Vec<Option<usize>> = vec![None; sequence.len()];
        let mut loop_neighbors = IntMap::default();
        let mut pair_list: IntMap<NAIDX, NAIDX>  = IntMap::default();
        let mut energy = 0;

        // Decomposing the structure into loops and initializing
        // loop_list, pair_list, and loop_lookup. 
        pairings.for_each_loop(|l| {
            let (lli, en) = registry.insert_loop(l.to_owned());
            if let Some((i, j)) = l.closing() {
                pair_list.insert(i as NAIDX, j as NAIDX); 
            }
            for k in &l.inclusive_unpaired_indices() {
                loop_lookup[*k] = Some(lli);
            }
            energy += en;
            let neighbors = registry.get_loop_neighbors(lli);
            loop_neighbors.insert(lli, neighbors);
        });

        let mut pair_neighbors = IntMap::default();
        for (i, j) in pair_list.iter() {
            let o_index = loop_lookup[*i as usize].unwrap();
            let i_index = loop_lookup[*j as usize].unwrap();
            // How does the free energy change if the move is applied.
            let delta = registry.calc_pair_energy(o_index, i_index);
            pair_neighbors.insert(*i, delta);
        }

        Ok(LoopStructure {
            registry,
            loop_lookup,
            loop_neighbors,
            pair_list,
            pair_neighbors,
            energy,
        })
    }

}

impl<'a, E: EnergyModel> From<&LoopStructure<'a, E>> for DotBracketVec {
    fn from(ls: &LoopStructure<'a, E>) -> Self {
        // Use the same logic as your Display impl, but avoid allocating a String unnecessarily
        let mut vec = vec![DotBracket::Unpaired; ls.registry.sequence.len()];
        for (i, j) in &ls.pair_list {
            vec[*i as usize] = DotBracket::Open;
            vec[*j as usize] = DotBracket::Close;
        }
        DotBracketVec(vec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;

    #[test]
    fn test_add_then_del_roundtrip() {
        let seq = NucleotideVec::from_lossy("GCCCCGGUCA");
        let structure = PairTable::try_from("..........").unwrap();
        let model = ViennaRNA::default();

        let mut ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();

        // Clone neighbor list so we don’t mutate while iterating
        let neighbors: Vec<(NAIDX, NAIDX, i32)> = ls
            .get_add_neighbors_per_loop()
            .iter()
            .flat_map(|(_, nbrs)| nbrs.iter().copied())
            .collect();

        for (i, j, de) in neighbors {
            let initial_energy = ls.energy();
            println!("({i} {j} {de}) at energy: {}", initial_energy);

            // add pair
            let _ = ls.apply_add_move(i, j);
            println!("{i} {j} {}", ls.energy());

            // delete the same pair
            let (p, q, rde) = ls.get_del_neighbors().first().cloned().unwrap();
            println!("({p} {q} {rde}) at energy: {}", ls.energy());
            assert_eq!((i, j), (p, q), "same pair gets deleted");
            assert_eq!(de, -rde, "inverse energy of reverse move");

            // delete pair 
            let _ = ls.apply_del_move(i, j);
            let roundtrip_energy = ls.energy();
            println!("{i} {j} {}", ls.energy());
            assert_eq!(roundtrip_energy, initial_energy, "roundtrip energy mismatch");
        }
    }

    #[test]
    fn test_add_then_del_bug() {
        let seq = NucleotideVec::from_lossy("GCCCCGGUCA");
        let structure = PairTable::try_from("((....).).").unwrap();
        let model = ViennaRNA::default();

        let ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();
        let neighbors = ls.get_del_neighbors();
        println!("{:?}", neighbors);

        let structure = PairTable::try_from("..........").unwrap();
        let mut ls = LoopStructure::try_from((&seq[..], &structure, &model)).unwrap();
        let _ = ls.apply_add_move(0, 8);
        println!("{:?}", neighbors);
        let _ = ls.apply_add_move(1, 6);
        assert_eq!(neighbors, ls.get_del_neighbors());
    }

    #[test]
    fn test_cotranscr() {
        let model = ViennaRNA::default();

        let seq = NucleotideVec::from_lossy("GCCCCGGUCA");
        let st1 = PairTable::try_from("(...)").unwrap();
        let st2 = PairTable::try_from("(...).").unwrap();

        let mut ls1 = LoopStructure::try_from((&seq[..], &st1, &model)).unwrap();
        assert_eq!(ls1.energy(), 540);
        let ls2 = LoopStructure::try_from((&seq[..], &st2, &model)).unwrap();
        assert_eq!(ls2.energy(), 370);

        ls1.apply_ext_move();
        assert_eq!(ls1.energy(), ls2.energy());
    }

    #[test]
    fn test_appyl_ext_move() {
        let seq = NucleotideVec::from_lossy("UGCCCCGGUCA");
        let structure_1 = PairTable::try_from(".((....).)").unwrap();
        let structure_2 = PairTable::try_from(".((....).).").unwrap();
        let model = ViennaRNA::default();

        let mut ls = LoopStructure::try_from((&seq[..], &structure_1, &model)).unwrap();
        let add_neighbors_1: Vec<_> = ls.get_add_neighbors_per_loop()
            .iter()
            .flat_map(|(_, neighbors) | neighbors.clone())
            .collect();

        let energy_1 = ls.energy();
        println!("Energy before extension: {}", ls.energy());
        

        let result = ls.apply_ext_move();

        let energy_2 = ls.energy();
        let add_neighbors_2: Vec<_> = ls.get_add_neighbors_per_loop()
            .iter()
            .flat_map(|(_, neighbors) | neighbors.clone())
            .collect();
        println!("Energy after extension:{}", ls.energy());
        assert_ne!(energy_1, energy_2);
        assert_ne!(add_neighbors_1, add_neighbors_2);

        let ls_exp = LoopStructure::try_from((&seq[..], &structure_2, &model)).unwrap();
        let energy_exp = ls_exp.energy();
        let add_neighbors_exp: Vec<_> = ls_exp.get_add_neighbors_per_loop()
            .iter()
            .flat_map(|(_, neighbors) | neighbors.clone())
            .collect();
        println!("Expected energy after extension:{}", energy_exp);
        assert_eq!(energy_2, energy_exp); //assert if energy after extension matches expected energy 
        assert_eq!(add_neighbors_2, add_neighbors_exp); //assert if neighbors match exptected neighbors
    }
    

}

