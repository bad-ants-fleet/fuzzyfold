use std::fmt;
use nohash_hasher::IntMap;
use nohash_hasher::IntSet;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::NearestNeighborLoop;
use ff_energy::LoopDecomposition;
use ff_energy::EnergyModel;
use ff_energy::Base;

struct LoopCache<'a, E: EnergyModel> {
    sequence: &'a [Base],
    model: &'a E,
    loop_list: IntMap<usize, (NearestNeighborLoop, i32)>,
    l_indices: IntSet<usize>,
}

impl<'a, E: EnergyModel> LoopCache<'a, E> {

    pub fn new(sequence: &'a [Base], model: &'a E) -> Self {
        Self { 
            sequence,
            model,
            loop_list: IntMap::default(),
            l_indices: IntSet::default(),
        }
    }

    pub fn insert_loop(&mut self, combo: &NearestNeighborLoop) -> usize {
        let energy = self.model.energy_of_loop(self.sequence, combo);
        let index = self.allocate_index();
        self.loop_list.insert(index, (combo.to_owned(), energy));
        index
    }

    pub fn get_loop_by_index(&self, index: &usize) -> &(NearestNeighborLoop, i32) {
        self.loop_list.get(index).expect("Incorrect index")
    }

    pub fn calc_pair_energy(&mut self, outer_index: usize, inner_index: usize) -> i32 {
        let (outer, o_en) = self.loop_list.get(&outer_index).expect("Missing outer loop_list entry.");
        let (inner, i_en) = self.loop_list.get(&inner_index).expect("Missing inner loop_list entry.");
        let combo = &outer.join_loop(inner);
        let combo_energy = self.model.energy_of_loop(self.sequence, combo);
        combo_energy - (o_en + i_en)
    }

    pub fn apply_delete_move(&mut self, outer_index: usize, inner_index: usize, delta: i32) -> usize {
        let (outer, o_en) = self.loop_list.get(&outer_index).expect("Missing outer loop_list entry.");
        let (inner, i_en) = self.loop_list.get(&inner_index).expect("Missing inner loop_list entry.");
        let combo = outer.join_loop(inner);
        let combo_energy = (o_en + i_en) - delta;

        // re-use outer_index for the new loop.
        self.loop_list.insert(outer_index, (combo, combo_energy));
        self.loop_list.remove(&inner_index);
        self.l_indices.insert(inner_index);
        outer_index
    }

    pub fn apply_addition_move(&mut self, 
        combo_index: usize, 
        combo: NearestNeighborLoop, c_energy: i32, i: NAIDX, j: NAIDX
    ) -> (usize, usize, i32) {
        let (outer, inner) = combo.split_loop(i, j);

        //NOTE: could look delta up directly by searching loop_list.
        let outer_energy = self.model.energy_of_loop(self.sequence, &outer);
        let inner_energy = self.model.energy_of_loop(self.sequence, &inner);

        let outer_index = combo_index;
        let inner_index = self.allocate_index();
        self.loop_list.insert(outer_index, (outer, outer_energy));
        self.loop_list.insert(inner_index, (inner, inner_energy));

        // How does the energy change if we apply the base-pair move.
        let delta = (outer_energy + inner_energy) - c_energy;
        (outer_index, inner_index, -delta)
    }

    fn get_loop_neighbors(&self, index: usize) -> MoveEnergies {
        let (combo, energy) = self.loop_list.get(&index).expect("where's the loop?");
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

    pub fn allocate_index(&mut self) -> usize {
        if let Some(&x) = self.l_indices.iter().next() {
            self.l_indices.remove(&x);
            x
        } else {
            self.loop_list.len()
        }
    }
}

type MoveEnergies = Vec<(NAIDX, NAIDX, i32)>;
type IndexedLoopNeighbors = (usize, MoveEnergies);

pub struct LoopStructure<'a, E: EnergyModel> {
    registry: LoopCache<'a, E>,
    /// From sequence index to registry index.
    loop_lookup: IntMap<NAIDX, usize>, 
    /// registry index to list of (i, j, deltaE)
    loop_neighbors: IntMap<usize, MoveEnergies>,
    /// Current pairs, i<j where i is the id.
    pair_list: IntMap<NAIDX, NAIDX>,
    /// pair id to deltaE
    pair_neighbors: IntMap<NAIDX, i32>, 
}

impl<'a, E: EnergyModel> Clone for LoopStructure<'a, E> {
    fn clone(&self) -> Self {
        Self {
            registry: LoopCache {
                sequence: self.registry.sequence,
                model: self.registry.model,
                loop_list: self.registry.loop_list.clone(),
                l_indices: self.registry.l_indices.clone(),
            },
            loop_lookup: self.loop_lookup.clone(),
            loop_neighbors: self.loop_neighbors.clone(),
            pair_list: self.pair_list.clone(),
            pair_neighbors: self.pair_neighbors.clone(),
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
    pub fn loop_lookup(&self) -> &IntMap<NAIDX, usize> {
        &self.loop_lookup
    }

    pub fn energy(&self) -> i32 {
        self.registry.loop_list
            .values()
            .map(|(_, e)| *e)
            .sum()
    }

    fn update_pair_neighbors(&mut self,
        pairs: &[(NAIDX, NAIDX)]
    ) -> MoveEnergies
    {
        let mut change = Vec::new();
        for &(i, j) in pairs {
            let &o_index = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");
            let &i_index = self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j.");
            let delta = self.registry.calc_pair_energy(o_index, i_index);
            self.pair_neighbors.insert(i, delta);
            change.push((i, j, delta));
        }
        change
    }
  
    pub fn apply_del_move(&mut self, i: NAIDX, j: NAIDX) -> 
        (MoveEnergies, MoveEnergies) {
        debug_assert_eq!(&j,
            self.pair_list.get(&i).expect("Missing pair_list entry."));
        let &delta = self.pair_neighbors.get(&i).expect("Missing pair_neighbors entry.");
        self.pair_list.remove(&i); 
        self.pair_neighbors.remove(&i); 

        let &o_index = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");
        let &i_index = self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j.");
        let c_id = self.registry.apply_delete_move(o_index, i_index, -delta);
        debug_assert_eq!(c_id, o_index); // by convention.

        let loop_neighbors = self.registry.get_loop_neighbors(o_index);
        self.loop_neighbors.insert(o_index, loop_neighbors.clone());
        self.loop_neighbors.remove(&i_index).expect("at least empty list.");

        let (combo, _) = self.registry.get_loop_by_index(&o_index);
        for k in &combo.inclusive_unpaired_indices() {
            debug_assert!(self.loop_lookup[&(*k as NAIDX)] == o_index || self.loop_lookup[&(*k as NAIDX)] == i_index);
            self.loop_lookup.insert(*k as NAIDX, o_index);
        }

        let pair_changes = self.update_pair_neighbors(&combo.pairs());
        (loop_neighbors, pair_changes)
    }

    pub fn apply_add_move(&mut self, i: NAIDX, j: NAIDX
    ) -> (
        MoveEnergies,
        MoveEnergies,
        MoveEnergies,
    ) {
        let &c_index = self.loop_lookup.get(&i).expect("Missing loop_lookup entry for i.");
        debug_assert_eq!(&c_index, 
            self.loop_lookup.get(&j).expect("Missing loop_lookup entry for j."),
            "Missing loop_lookup entry for j."
        );
        let (combo, c_energy) = self.registry.get_loop_by_index(&c_index);
        let combo_pairs = &combo.pairs();
        // How does the energy change if we apply the base-pair move.
        let (o_id, i_id, delta) = self.registry.apply_addition_move(
            c_index, combo.clone(), *c_energy, i, j
        );


        let new_outer_add_neighbors = self.registry.get_loop_neighbors(o_id);
        let new_inner_add_neighbors = self.registry.get_loop_neighbors(i_id);

        self.loop_neighbors.insert(o_id, new_outer_add_neighbors.clone());
        self.loop_neighbors.insert(i_id, new_inner_add_neighbors.clone());
        self.pair_list.insert(i, j);
        self.pair_neighbors.insert(i, delta);

        let (outer, _) = self.registry.get_loop_by_index(&o_id);
        for k in &outer.inclusive_unpaired_indices() {
            self.loop_lookup.insert(*k as NAIDX, o_id);
        }
        let (inner, _) = self.registry.get_loop_by_index(&i_id);
        for k in &inner.inclusive_unpaired_indices() {
            self.loop_lookup.insert(*k as NAIDX, i_id);
        }

        let mut pair_changes = self.update_pair_neighbors(combo_pairs);
        pair_changes.push((i, j, delta));

        (new_outer_add_neighbors, new_inner_add_neighbors, pair_changes)
    }

}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a [Base], &T, &'a E)> for LoopStructure<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a [Base], &T, &'a E)
    ) -> Result<Self, Self::Error> {
        let mut registry = LoopCache::new(sequence, model);
        let mut pair_list: IntMap<NAIDX, NAIDX>  = IntMap::default();
        let mut loop_lookup: IntMap<NAIDX, usize> = IntMap::default();

        // Decomposing the structure into loops and initializing
        // loop_list, pair_list, and loop_lookup. 
        pairings.for_each_loop(|l| {
            let lli = registry.insert_loop(l);
            if let Some((i, j)) = l.closing() {
                pair_list.insert(i as NAIDX, j as NAIDX); 
            }
            for k in &l.inclusive_unpaired_indices() {
                loop_lookup.insert(*k as NAIDX, lli);
            }
        });

        // Now we want to get all neighbors where pairs can be added
        let mut loop_neighbors = IntMap::default();
        for &nn_idx in registry.loop_list.keys() {
            let neighbors = registry.get_loop_neighbors(nn_idx);
            loop_neighbors.insert(nn_idx, neighbors);
        }

        let mut pair_neighbors = IntMap::default();
        for (i, j) in pair_list.iter() {
            let &o_index = loop_lookup.get(i).expect("Missing loop_lookup entry for i.");
            let &i_index = loop_lookup.get(j).expect("Missing loop_lookup entry for j.");
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

}

