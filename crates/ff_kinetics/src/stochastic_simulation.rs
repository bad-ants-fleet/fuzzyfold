use std::fmt;
use rand::Rng; // -> R
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_energy::EnergyModel; // -> E

use crate::RateModel; // -> K
use crate::LoopStructure;
use crate::explore::Move;
use crate::rate_tree::RateTree;
use crate::rate_tree::RateList;

/// An SSA implementation for LoopStructure.
pub struct LoopStructureSSA<'a, E: EnergyModel, K: RateModel> {
    /// The current RNA structure representation.
    loopstructure: LoopStructure<'a, E>,
    /// Anything with the RateModel trait.
    ratemodel: &'a K,
    /// Heap-like data structure for sampling.
    rate_tree: RateTree,
    /// Selecting the move from loop-reaction.
    loop_lists: IntMap<usize, RateList>,
}

impl<'a, E, K> fmt::Debug for LoopStructureSSA<'a, E, K>
where
    E: EnergyModel,
    K: RateModel + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructureSSA")
            .field("ratemodel", &self.ratemodel) 
            .field("loopstructure", &format!("{}", self.loopstructure))
            .finish()
    }
}

impl<'a, E: EnergyModel, K: RateModel> From<(LoopStructure<'a, E>, &'a K)>
    for LoopStructureSSA<'a, E, K>
{
    fn from((loopstructure, ratemodel): (LoopStructure<'a, E>, &'a K)) -> Self {
        let mut rate_tree = RateTree::new(loopstructure.len());
        let mut loop_lists = IntMap::default();

        for (lli, add_neighbors) in loopstructure.get_add_neighbors_per_loop().iter() {
            let mut llist = RateList::new(add_neighbors.len());
            for &(i, j, delta) in add_neighbors {
                llist.insert(Move::Add { i, j }, ratemodel.rate(delta));
            }
            if !llist.is_empty() {
                assert!(*lli < 65000, "can't convert to NAIDX!");
                rate_tree.insert(Move::Loop { idx: *lli as NAIDX }, llist.total_rate())
            }
            loop_lists.insert(*lli, llist);
        }

        for (i, j, delta) in loopstructure.get_del_neighbors() {
            rate_tree.insert(Move::Del { i, j }, ratemodel.rate(delta));
        }

        Self {
            ratemodel,
            loopstructure,
            rate_tree,
            loop_lists,
        }
    }
}

impl<'a, E: EnergyModel, K: RateModel> LoopStructureSSA<'a, E, K> {
    pub fn current_structure(&self) -> String {
        format!("{}", self.loopstructure)
    }   

    pub fn remove_loop_reaction(&mut self, i: NAIDX) {
        let lli = self.loopstructure.loop_lookup().get(&i).unwrap();
        let llist = self.loop_lists.remove(lli).expect("Reaction must exist.");
        if !llist.is_empty() && !self.rate_tree.remove(&Move::Loop { idx: *lli as NAIDX}) {
            panic!("Could not find loop-reaction in RateTree!");
        }
    }

    pub fn remove_pair_reaction(&mut self, (i, j): (NAIDX, NAIDX)) {
        if !self.rate_tree.remove(&Move::Del { i, j }) {
            panic!("Could not find pair-reaction in RateTree!");
        }
        self.remove_loop_reaction(i);
        self.remove_loop_reaction(j);
    }

    pub fn insert_loop_reactions(&mut self, 
        lli: usize, 
        add_neighbors: Vec<(NAIDX, NAIDX, i32)>
    ) {
        let mut llist = RateList::new(add_neighbors.len());
        for (i, j, delta) in add_neighbors {
            llist.insert(Move::Add { i, j }, self.ratemodel.rate(delta));
        }
        if !llist.is_empty() {
            assert!(lli < 65000, "can't convert to NAIDX!");
            self.rate_tree.insert(Move::Loop { idx: lli as NAIDX }, llist.total_rate())
        }
        self.loop_lists.insert(lli, llist);
    }

    pub fn update_pair_reactions(&mut self, change: Vec<(NAIDX, NAIDX, i32)>) {
        for (i, j, delta) in change {
            let mv = Move::Del { i, j };
            let rate = self.ratemodel.rate(delta);
            if !self.rate_tree.update_rate(&mv, rate) {
                self.rate_tree.insert(mv, rate);
            }
        }
    }

    pub fn simulate<R, F>(
        &mut self,
        rng: &mut R,
        t_max: f64,
        mut callback: F,
    )
    where
        R: Rng + ?Sized,
        F: FnMut(f64, f64, f64, &LoopStructure<'a, E>) -> bool,
    {
        let mut t = 0.;

        while t < t_max {
            let rsum = self.rate_tree.total_rate();

            // sample waiting time ~ Exp(flux)
            let tinc = -rng.random::<f64>().ln() / rsum;

            // Callback bewore applying the waiting time.
            // If callback return's false, then abort the simulation!
            if !callback(t, tinc, rsum, &self.loopstructure) {
                break;
            }

            t += tinc;

            let threshold = rng.random::<f64>() * rsum;
            let mv = self.rate_tree.select_by_threshold(threshold);

            match mv {
                Some((_, Move::Del { i, j })) => {
                    self.remove_pair_reaction((i, j));
                    let ((lli, neighbors), pair_changes) = self
                        .loopstructure.apply_del_move(i, j);
                    self.insert_loop_reactions(lli, neighbors);
                    self.update_pair_reactions(pair_changes);
                },
                Some((th, Move::Loop { idx })) => {
                    let llist = self.loop_lists.get(&(idx as usize)).unwrap();
                    let (i, j) = llist.select_by_threshold(th).ij();
                    self.remove_loop_reaction(i);
                    let ((lli, ami), (llj, amj), pair_changes) = self
                        .loopstructure.apply_add_move(i, j);
                    self.insert_loop_reactions(lli, ami);
                    self.insert_loop_reactions(llj, amj);
                    self.update_pair_reactions(pair_changes);
                },
                _ => panic!("No reaction chosen despite positive flux"),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;

    // --- The actual test ----------------------------------------------------
    #[test]
    fn test_simple_ssa_simulation() {
        let emodel = ViennaRNA::default();
        let rmodel = Metropolis::new(emodel.temperature(), 1.0);
        let mut rng = StdRng::seed_from_u64(42);

        let sequence = "CAAAG";
        let structure = ".....";

        let sequence = NucleotideVec::try_from(sequence).unwrap();
        let pairings = PairTable::try_from(structure)
            .expect("invalid structure in input");
        let loops = LoopStructure::try_from((&sequence[..], &pairings, &emodel))
            .expect("failed to build loop structure");

        let mut simulator = LoopStructureSSA::from((loops, &rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng, 1.0, |t, tinc, flux, _ls| {
            steps += 1;
            println!(
                "Step {}: t={:.4}, Δt={:.4}, flux={:.3e}",
                steps, t, tinc, flux
            );
            true
        });

        assert!(steps > 0, "Simulation must perform at least one step");
        //assert!(simulator.log_flux.is_finite(), "Flux must remain finite");
    }
}



