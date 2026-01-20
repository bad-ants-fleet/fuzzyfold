use ff_structure::DotBracketVec;
use rand::Rng; // -> R

use crate::Walker;
use crate::RateModel; // -> K
use crate::rate_tree::RateTree;

/// An SSA implementation for LoopStructure.
pub struct SSA<'a, W: Walker, K: RateModel> {
    /// The current RNA structure representation.
    walker: W,
    /// Anything with the RateModel trait.
    ratemodel: &'a K,
    /// Heap-like data structure for sampling.
    rate_tree: RateTree,
}

impl<'a, W: Walker, K: RateModel> 
From<(W, &'a K)> for SSA<'a, W, K>
{
    fn from((walker, ratemodel): (W, &'a K)) -> Self {
        let mut rate_tree = RateTree::new(walker.len());

        for (mv, delta) in walker.propose_moves() {
            let k = ratemodel.rate(&mv, delta);
            rate_tree.init_insert(mv, k);
        }
        rate_tree.init_partial_sums();

        Self {
            walker,
            ratemodel,
            rate_tree,
        }
    }
}

impl<'a, W: Walker, K: RateModel> SSA<'a, W, K> {
    pub fn current_structure(&self) -> DotBracketVec {
        self.walker.current_structure()
    }   

    pub fn current_energy(&self) -> i32 {
        self.walker.current_energy()
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
        F: FnMut(f64, f64, f64, &W) -> bool,
    {
        let mut t = 0.;

        while t < t_max {
            let rsum = self.rate_tree.total_rate();

            // sample waiting time ~ Exp(flux)
            let tinc = -rng.random::<f64>().ln() / rsum;

            // Callback bewore applying the waiting time.
            // If callback return's false, then abort the simulation!
            if !callback(t, tinc, rsum, &self.walker) {
                break;
            }

            t += tinc;

            let threshold = rng.random::<f64>() * rsum;
            let mv = self.rate_tree.select_by_threshold(threshold).expect("Must select a move!");
            let (old, new) = self.walker.apply_move(&mv);

            //println!("from: {:?} to: {:?} ", mv, mv.inverse());
            if !self.rate_tree.dirty_replace(&mv, &mv.inverse()) {
                panic!("Dirty rate replacement failed.");
            }

            for (mv, _) in old {
                //println!("remove: {:?}", mv);
                self.rate_tree.remove(mv);
            }

            for (mv, delta) in new {
                let k = self.ratemodel.rate(&mv, delta);
                if !self.rate_tree.update_rate(&mv, k) {
                    //println!("insert: {:?} -> {}", mv, k);
                    self.rate_tree.insert(mv, k);
                } 
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
    use ff_energy::EnergyModel;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;
    use crate::movesets::AddDelMoves;
    use crate::movesets::LoopTable;

    #[test]
    fn test_simple_ssa_simulation() {
        let emodel = ViennaRNA::default();
        let rmodel = Metropolis::new(emodel.temperature(), 1.0);
        let mut rng = StdRng::seed_from_u64(42);

        let sequence = "CAAAG";
        let pairings = ".....";

        let sequence = NucleotideVec::try_from(sequence).unwrap();
        let pairings = PairTable::try_from(pairings).unwrap();

        let ltab1 = LoopTable::try_from((&sequence[..], &pairings, &emodel)).unwrap();
        let adm1 = AddDelMoves::from(ltab1);

        let mut simulator = SSA::from((adm1, &rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng, 1.0, |t, tinc, flux, _| {
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



