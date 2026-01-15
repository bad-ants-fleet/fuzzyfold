use rand::Rng; // -> R

use ff_structure::NAIDX;
use ff_energy::EnergyModel; // -> E

use crate::RateModel; // -> K
use crate::LoopStructure;
use crate::explore::Move;
use crate::rate_tree::RateTree;

/// An SSA implementation for LoopStructure.
pub struct LoopStructureSSA<'a, E: EnergyModel, K: RateModel> {
    /// The current RNA structure representation.
    loopstructure: LoopStructure<'a, E>,
    /// Anything with the RateModel trait.
    ratemodel: &'a K,
    /// Heap-like data structure for sampling.
    rate_tree: RateTree,
}

// NOTE: A basic implementation for debugging. Feel free to adapt.
// This code is not performance-critical and not considered part of the stable API.
#[cfg(debug_assertions)]
impl<'a, E, K> std::fmt::Debug for LoopStructureSSA<'a, E, K>
where
    E: EnergyModel,
    K: RateModel + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LoopStructureSSA")
            .field("loopstructure", &format!("{}", self.loopstructure))
            .field("ratemodel", &self.ratemodel) 
            .finish()
    }
}

impl<'a, E: EnergyModel, K: RateModel> From<(LoopStructure<'a, E>, &'a K)>
    for LoopStructureSSA<'a, E, K>
{
    fn from((loopstructure, ratemodel): (LoopStructure<'a, E>, &'a K)) -> Self {
        let mut rate_tree = RateTree::new(loopstructure.len());

        for (i, j, delta) in loopstructure.iter_add_moves() {
            rate_tree.init_insert(
                Move::Add { i, j }, 
                ratemodel.rate(delta)
            );
        }

        for (i, j, delta) in loopstructure.iter_del_moves() {
            rate_tree.init_insert(
                Move::Del { i, j },
                ratemodel.rate(delta)
            );
        }

        rate_tree.init_partial_sums();

        Self {
            ratemodel,
            loopstructure,
            rate_tree,
        }
    }
}

impl<'a, E: EnergyModel, K: RateModel> LoopStructureSSA<'a, E, K> {
    pub fn current_structure(&self) -> &LoopStructure<'a, E> {
        &self.loopstructure
    }   

    /// Main simulation function.
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
                Some(del @ Move::Del { i, j }) => {
                    // We use "dirty replace, without changing the rate, because the loop 
                    // will be updated afterwards as part of the new loop reactions.
                    if !self.rate_tree.dirty_replace(&del, &Move::Add {i, j}) {
                        panic!("Dirty rate replacement failed.");
                    }
                    let (neighbors, pair_changes) = self
                        .loopstructure.apply_del_move(i, j);
                    self.update_loop_reactions(neighbors);
                    self.update_pair_reactions(pair_changes);
                },
                Some(add @ Move::Add { i, j }) => {
                    // We use "dirty replace, without changing the rate, because the loop 
                    // will be updated afterwards as part of the new pair reactions.
                    if !self.rate_tree.dirty_replace(&add, &Move::Del {i, j}) {
                        panic!("Dirty rate replacement failed.");
                    }
                    // Get the loop-list index to remove loop reactions.
                    let lli = self.loopstructure.loop_lookup()[i as usize];
                    for (p, q, _) in self.loopstructure.iter_add_moves_for_loop(lli) {
                        // Those are the ones that will be updated later anyway.
                        if q < i || j < p || (i < p && q < j) || (p < i && j < q) {
                            continue;
                        }
                        self.rate_tree.remove(Move::Add { i: p, j: q });
                    }
                    let (ami, amj, pair_changes) = self
                        .loopstructure.apply_add_move(i, j);
                    self.update_loop_reactions(ami);
                    self.update_loop_reactions(amj);
                    self.update_pair_reactions(pair_changes);
                },
                None => panic!("No reaction chosen despite positive flux"),
            }
        }
    }

    /// Update all reactions to add a new pair within a specific loop. 
    fn update_loop_reactions(
        &mut self, 
        add_neighbors: Vec<(NAIDX, NAIDX, i32)>
    ) {
        for (i, j, delta) in add_neighbors {
            let mv = Move::Add { i, j };
            let rate = self.ratemodel.rate(delta);
            if !self.rate_tree.update_rate(&mv, rate) {
                self.rate_tree.insert(mv, rate);
            }
        }
    }

    /// Update all reactions that remove an existing pair. 
    fn update_pair_reactions(
        &mut self, 
        change: Vec<(NAIDX, NAIDX, i32)>
    ) {
        for (i, j, delta) in change {
            let mv = Move::Del { i, j };
            let rate = self.ratemodel.rate(delta);
            if !self.rate_tree.update_rate(&mv, rate) {
                self.rate_tree.insert(mv, rate);
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

    #[test]
    fn test_simple_ssa_simulation() {
        let emodel = ViennaRNA::default();
        let rmodel = Metropolis::new(emodel.temperature(), 1.0);
        let mut rng = StdRng::seed_from_u64(42);

        let sequence = "CAAAG";
        let pairings = ".....";

        let sequence = NucleotideVec::try_from(sequence).unwrap();
        let pairings = PairTable::try_from(pairings).unwrap();
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
    }
}



