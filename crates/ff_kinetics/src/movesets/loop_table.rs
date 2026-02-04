use std::fmt;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_energy::EnergyModel;
use ff_energy::NucleotideVec;
use ff_energy::LoopDecomposition;
use ff_energy::NearestNeighborLoop;

type LoopEntry = (NearestNeighborLoop, i32);

/// Direct storage of the loop decomposition from a structure.
///
/// Keeps all loops in a vector, and can efficiently look up
/// the loops corresponding to paired and unpaired positions.
///
/// Keeps track of base-pairs and the current energy for convenience.
///
pub struct LoopTable<'a, E: EnergyModel> {
    sequence: &'a NucleotideVec,
    model: &'a E,
    loops: Vec<LoopEntry>,
    stale: Vec<usize>,
    loop_lookup: Vec<usize>,
    pair_lookup: IntMap<NAIDX, NAIDX>,
    energy: i32,
}

/// Clone creates an independent loop bookkeeping state,
/// but shares the underlying sequence and energy model.
impl<'a, E: EnergyModel> Clone for LoopTable<'a, E> {
    fn clone(&self) -> Self {
        Self {
            sequence: self.sequence,
            model: self.model,
            loops: self.loops.clone(),
            stale: self.stale.clone(),
            loop_lookup: self.loop_lookup.clone(),
            pair_lookup: self.pair_lookup.clone(),
            energy: self.energy,
        }
    }
}

impl<'a, E: EnergyModel> LoopTable<'a, E> {
    pub fn min_hairpin_size(&self) -> usize {
        self.model.min_hairpin_size()
    }

    pub fn can_pair(&self, i: usize, j: usize) -> bool {
        self.model.can_pair(self.sequence[i], self.sequence[j])
    }

    pub fn energy_of_loop(&self, nn_loop: &NearestNeighborLoop) -> i32 {
        self.model.energy_of_loop(self.sequence, nn_loop)
    }

    pub fn loops_len(&self) -> usize {
        self.loops.len()
    }

    pub fn lookup_len(&self) -> usize {
        self.loop_lookup.len()
    }

    pub fn pairs(&self) -> impl Iterator<Item = (&NAIDX, &NAIDX)> + '_ {
        self.pair_lookup.iter()
    }

    pub fn delete_pair(&mut self, i: &NAIDX) -> NAIDX {
        self.pair_lookup.remove(i).expect("Missing pair_list entry.")
    }

    pub fn insert_pair(&mut self, i: NAIDX, j: NAIDX) {
        self.pair_lookup.insert(i, j);
    }

    pub fn energy(&self) -> i32 {
        self.energy
    }

    pub fn update_lookup(&mut self, idx: usize) {
        let (nn_loop, _) = self.get(idx);
        for k in nn_loop.inclusive_unpaired_indices() {
            self.loop_lookup[k] = idx;
        }
    }

    pub fn pair_lookup(&self, idx: &NAIDX) -> NAIDX {
        self.pair_lookup[idx]
    }

    pub fn loop_lookup(&self, idx: usize) -> usize {
        self.loop_lookup[idx]
    }

    pub fn geti(&self, i: usize) -> &LoopEntry {
        &self.loops[self.loop_lookup[i]]
    }

    pub fn get(&self, idx: usize) -> &LoopEntry {
        self.loops
            .get(idx)
            .expect("Invalid LoopCache index")
    }

    pub fn mark_stale(&mut self, idx: usize) {
        self.energy -= self.loops[idx].1;
        self.stale.push(idx);
    }

    pub fn insert_loopentry(
        &mut self, 
        index: Option<usize>, 
        entry: LoopEntry,
    ) -> usize {
        self.energy += entry.1;
        if let Some(idx) = index { 
            self.energy -= self.loops[idx].1;
            self.loops[idx] = entry;
            idx
        } else if let Some(idx) = self.stale.pop() {
            self.loops[idx] = entry;
            idx
        } else {
            let idx = self.loops.len();
            self.loops.push(entry);
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

impl<'a, E: EnergyModel> From<&LoopTable<'a, E>> for DotBracketVec {
    fn from(ltab: &LoopTable<'a, E>) -> Self {
        // Use the same logic as your Display impl, but avoid allocating a String unnecessarily
        let mut vec = vec![DotBracket::Unpaired; ltab.lookup_len()];
        for (i, j) in ltab.pairs() {
            vec[*i as usize] = DotBracket::Open;
            vec[*j as usize] = DotBracket::Close;
        }
        DotBracketVec(vec)
    }
}

impl<'a, T: LoopDecomposition, E: EnergyModel> TryFrom<(&'a NucleotideVec, &T, &'a E)> 
for LoopTable<'a, E> {
    type Error = String;

    fn try_from((sequence, pairings, model): (&'a NucleotideVec, &T, &'a E)
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

