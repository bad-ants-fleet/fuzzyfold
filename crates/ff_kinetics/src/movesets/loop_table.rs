use std::fmt;
use std::sync::Arc;
use nohash_hasher::IntMap;

use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::MultiPairTable;
use ff_energy::Base;
use ff_energy::EnergyModel;
use ff_energy::Gap;
use ff_energy::NucleotideVec;
use ff_energy::NestedLoop;
use ff_energy::Ring;
use ff_energy::RingDecomposition;

use crate::MoveSetError;


pub(crate) type LoopEntry = (NestedLoop, i32);

/// Direct storage of the loop decomposition from a structure.
///
/// Keeps all loops in a vector, and can efficiently look up
/// the loops corresponding to paired and unpaired positions.
///
/// Keeps track of base-pairs and the current energy for convenience.
///
pub struct LoopTable<E: EnergyModel> {
    sequence: NucleotideVec,
    model: Arc<E>,
    loops: Vec<LoopEntry>,
    stale: Vec<usize>,
    loop_lookup: Vec<usize>,
    pair_lookup: IntMap<usize, usize>,
    energy: i32,
}

/// Clone creates an independent loop bookkeeping state,
/// but shares the underlying sequence and energy model.
impl<E: EnergyModel> Clone for LoopTable<E> {
    fn clone(&self) -> Self {
        Self {
            sequence: self.sequence.clone(),
            model: self.model.clone(),
            loops: self.loops.clone(),
            stale: self.stale.clone(),
            loop_lookup: self.loop_lookup.clone(),
            pair_lookup: self.pair_lookup.clone(),
            energy: self.energy,
        }
    }
}

impl<E: EnergyModel> LoopTable<E> {
    pub fn sequence_length(&self) -> usize {
        self.sequence.len()
    }

    pub fn min_hairpin_size(&self) -> usize {
        self.model.min_hairpin_size()
    }

    /// True for a circular strand: its last position is the strand wrap.
    pub fn is_circular(&self) -> bool {
        self.sequence.last() == Some(&Base::SW)
    }

    pub fn can_pair(&self, i: usize, j: usize) -> bool {
        self.model.can_pair(self.sequence[i], self.sequence[j])
    }

    /// False for a hairpin with fewer unpaired bases than the model allows.
    /// Index distances rule these out, except for the hairpin through the 
    /// strand wrap of a circular strand.
    pub(crate) fn fits(&self, l: &NestedLoop) -> bool {
        l.hairpin_size().is_none_or(|n| n >= self.min_hairpin_size())
    }

    pub(crate) fn energy_of_loop(&self, l: &NestedLoop) -> i32 {
        self.model.energy_of_nested_loop(&self.sequence, l)
            .expect("Broken energy evaluation!")
    }

    pub(crate) fn loops_len(&self) -> usize {
        self.loops.len()
    }

    pub fn lookup_len(&self) -> usize {
        self.loop_lookup.len()
    }

    pub fn pairs(&self) -> impl Iterator<Item = (&usize, &usize)> + '_ {
        self.pair_lookup.iter()
    }

    pub(crate) fn delete_pair(&mut self, i: &usize) -> usize {
        self.pair_lookup.remove(i).expect("Missing pair_list entry.")
    }

    pub(crate) fn insert_pair(&mut self, i: usize, j: usize) {
        self.pair_lookup.insert(i, j);
    }

    pub fn energy(&self) -> i32 {
        self.energy
    }

    pub(crate) fn extend_lookup(&mut self, idx: usize) {
        self.loop_lookup.push(idx);
    }

    pub(crate) fn update_lookup(&mut self, idx: usize) {
        let (l, _) = self.loops.get(idx).expect("Invalid LoopCache index");
        for s in l.stretches() {
            self.loop_lookup[s.unpaired.start..=s.ends.1].fill(idx);
        }
    }

    pub(crate) fn pair_lookup(&self, idx: &usize) -> usize {
        self.pair_lookup[idx]
    }

    pub(crate) fn loop_lookup(&self, idx: usize) -> usize {
        self.loop_lookup[idx]
    }

    pub(crate) fn geti(&self, i: usize) -> &LoopEntry {
        &self.loops[self.loop_lookup[i]]
    }

    pub(crate) fn get(&self, idx: usize) -> &LoopEntry {
        self.loops
            .get(idx)
            .expect("Invalid LoopCache index")
    }

    pub(crate) fn mark_stale(&mut self, idx: usize) {
        self.energy -= self.loops[idx].1;
        self.stale.push(idx);
    }

    pub(crate) fn insert_loopentry(
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

impl<E: EnergyModel> fmt::Display for LoopTable<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Convert sequence to string
        let mut dbr = vec!['.'; self.loop_lookup.len()];
        for (i, j) in &self.pair_lookup {
            dbr[*i ] = '(';
            dbr[*j ] = ')';
        }
        if self.is_circular() {
            *dbr.last_mut().unwrap() = '^';
        }
        let dbr_str: String = dbr.into_iter().collect();
        write!(f, "{}", dbr_str)
    }
}

impl<E: EnergyModel> fmt::Debug for LoopTable<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopTable")
            .field("loop_lookup", &self.loop_lookup)
            .field("pair_lookup", &self.pair_lookup)
            .finish()
    }
}

impl<E: EnergyModel> From<&LoopTable<E>> for DotBracketVec {
    fn from(ltab: &LoopTable<E>) -> Self {
        // Use the same logic as your Display impl, but avoid allocating a String unnecessarily
        let mut vec = vec![DotBracket::UNPAIRED; ltab.lookup_len()];
        for (i, j) in ltab.pairs() {
            vec[*i ] = DotBracket::OPEN;
            vec[*j ] = DotBracket::CLOSE;
        }
        if ltab.is_circular() {
            *vec.last_mut().unwrap() = DotBracket::STRAND_WRAP;
        }
        DotBracketVec(vec)
    }
}

impl<E: EnergyModel> TryFrom<(NucleotideVec, &MultiPairTable, Arc<E>)> for LoopTable<E> {
    type Error = MoveSetError;

    fn try_from((sequence, pairings, model): (NucleotideVec, &MultiPairTable, Arc<E>)
    ) -> Result<Self, Self::Error> {
        Self::new(sequence, pairings, model)
    }
}

impl<E: EnergyModel> LoopTable<E> {
    /// The table of a structure.
    pub fn new(sequence: NucleotideVec, pairings: &MultiPairTable, model: Arc<E>
    ) -> Result<Self, MoveSetError> {
        // One strand, linear or circular: a circular strand ends with the
        // strand wrap, the only strand end of the sequence.
        let circular = sequence.has_circular_strand();
        let strand_ends = sequence.iter().filter(|b| matches!(b, Base::SB | Base::SW)).count();
        if circular && (strand_ends > 1 || sequence.last() != Some(&Base::SW)) {
            return Err(MoveSetError::Unsupported);
        }

        let mut loops = Vec::new();
        let mut loop_lookup = Vec::with_capacity(sequence.len());
        let mut pair_lookup: IntMap<usize, usize>  = IntMap::default();
        let mut energy = 0;

        let mut current_len = 0;
        let mut exterior = 0;
        // Crossing pairs, also those that close no cycle, make a pseudoknot.
        let mut unsupported = !pairings.is_nested();
        let mut error = None;
        pairings.for_each_ring(|r| {
            // Kinetics support one strand without pseudoknots: one exterior
            // loop (through the wrap for a circular strand at full length),
            // and no joint exterior loops.
            let Ring::Loop(l) = r else {
                unsupported = true;
                return;
            };
            match l.gap() {
                Gap::Open { five, three } if !circular && five <= three => exterior += 1,
                Gap::Wrap { .. } => exterior += 1,
                Gap::None => (),
                _ => unsupported = true,
            }
            if unsupported || exterior > 1 || error.is_some() {
                return;
            }
            let (_, b) = l.span();
            let b = b + 1;
            if b > current_len {
                loop_lookup.resize(b, usize::MAX);
                current_len = b;
            }

            let loop_energy = match model.energy_of_nested_loop(&sequence, l) {
                Ok(en) => en,
                Err(e) => {
                    error = Some(e);
                    return;
                }
            };
            energy += loop_energy;

            if let Some((i, j)) = l.closing() {
                pair_lookup.insert(i, j); 
            }

            let loop_index = loops.len();
            for s in l.stretches() {
                loop_lookup[s.unpaired.start..=s.ends.1].fill(loop_index);
            }
            loops.push((l.to_owned(), loop_energy));
        });

        if let Some(e) = error {
            return Err(MoveSetError::Energy(e));
        }
        if unsupported || exterior != 1 {
            return Err(MoveSetError::Unsupported);
        }
        if circular {
            // The strand wrap belongs to no loop.
            loop_lookup.push(usize::MAX);
        }
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

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::Fuzzyfold;
    use ff_structure::MultiPairTable;

    #[test]
    fn test_one_strand_without_pseudoknots() {
        let model = Arc::new(Fuzzyfold::default());
        let table = |seq: &str, db: &str| LoopTable::try_from((
            NucleotideVec::try_from(seq).unwrap(),
            &MultiPairTable::try_from(db).unwrap(),
            model.clone(),
        ));
        assert!(table("GGGAAACCC", "(((...)))").is_ok());
        assert!(table("GGGAAACCCAAA", "((([[[)))]]]").is_err());
        assert!(table("GGG+CCC", "(((+)))").is_err());
        assert!(table("GGG+CCC", "...+...").is_err());
        // A circular strand alone, at full length.
        assert!(table("GGGGAAACCCCAAAAAA^", "((((...))))......^").is_ok());
        assert!(table("GGGGAAACCCCAAAAAA^", "((((...))))......").is_err());
        assert!(table("GGGGAAACCCCAAAAAA^AA", "((((...))))......^..").is_err());
        // The hairpin through the wrap has no unpaired bases.
        assert!(table("GGGAAACCC^", "(((...)))^").is_err());
    }
}
