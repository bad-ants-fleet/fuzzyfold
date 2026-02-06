use ff_structure::DotBracketVec;
use rand::Rng; // -> R

use crate::Walker;
use crate::Moves;
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
        let mut rate_tree = RateTree::default();

        for (mv, delta) in walker.propose_moves() {
            let k = ratemodel.rate(&mv, delta);
            if k > 0.0 {
                rate_tree.init_insert(mv, k);
            }
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

    pub fn co_simulate<R, F>(
        &mut self,
        rng: &mut R,
        times: &[f64],
        mut callback: F,
    ) where
        R: Rng + ?Sized,
        F: FnMut(f64, f64, f64, &W) -> bool,
    {

        //NOTE: This panic is temporary, should return an error eventually.
        if self.walker.sequence_length() != 
            self.walker.current_structure_length() - 1 + times.len() {
            panic!("mismatch between simulation times and sequence length");
        }

        let mut gtime = 0.0;
        for (idx, &time) in times.iter().enumerate() {
            // Wrap the user callback
            let mut co_callback = |t: f64, tinc: f64, 
                rsum: f64, w: &W| {
                    callback(t + gtime, tinc, rsum, w)
            };
            self.simulate(rng, time, &mut co_callback);

            // Skip extension after the last time point
            if idx + 1 == times.len() {
                break;
            }

            let (old, new) = self.walker.apply_extension(); 
            self.update_rate_tree(old, new);

            gtime += time;
        }
    }

    //NOTE: can be useful for debugging.
    fn _simple_update_rate_tree(&mut self, old: Moves, new: Moves) {
        for (mv, _) in old {
            self.rate_tree.remove(mv);

        }
        for (mv, delta) in new {
            let k = self.ratemodel.rate(&mv, delta);
            if k > 0. && !self.rate_tree.update_rate(&mv, k) {
                self.rate_tree.insert(mv, k);
            } 
        }
    }

    fn update_rate_tree(&mut self, old: Moves, new: Moves) {
        let mut del = old.iter();
        let mut add = new.iter();
        let mut cur_del = del.next();
        let mut cur_add = add.next();
        while cur_del.is_some() || cur_add.is_some() {
            match (cur_del, cur_add) {
                (Some((omv, _)), Some((nmv, delta))) => {
                    let k = self.ratemodel.rate(nmv, *delta);
                    if k == 0.0 || self.rate_tree.update_rate(nmv, k) {
                        cur_add = add.next();
                    } else if self.rate_tree.replace(omv, nmv, k) {
                        // Only true if the old move exists!
                        // This allows us to exlude valid moves with rate 0!
                        cur_add = add.next();
                        cur_del = del.next();
                    } else {
                        cur_del = del.next();
                    }
                }
                (Some((omv, _)), None) => {
                    self.rate_tree.remove(*omv);
                    cur_del = del.next();
                }
                (None, Some((nmv, delta))) => {
                    let k = self.ratemodel.rate(nmv, *delta);
                    if k > 0. && !self.rate_tree.update_rate(nmv, k) {
                        self.rate_tree.insert(*nmv, k);
                    }
                    cur_add = add.next();
                }
                (None, None) => unreachable!(),
            }
        }
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

            if rsum == 0.0 {
                callback(t, t_max, rsum, &self.walker);
                break;
            }

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
            self.update_rate_tree(old, new);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use rand::rng;
    use std::collections::HashSet;

    use ff_structure::PairTable;
    use ff_energy::ViennaRNA;
    use ff_energy::EnergyModel;
    use ff_energy::NucleotideVec;
    use crate::Metropolis;
    use crate::movesets::LoopNeighbors;
    use crate::movesets::shift_policy;
    use crate::movesets::loop_table::LoopTable;
    use crate::movesets::Move;

    macro_rules! setup_ssa_input {
        ($wname:ident, $rname:ident, $seq:expr, $db:expr) => {
            let emodel = ViennaRNA::default();
            let $rname = Metropolis::new(emodel.temperature(), 1.0, None, None);

            let sequence = NucleotideVec::try_from($seq)
                .expect("Invalid sequence?");
            let pairings = PairTable::try_from($db)
                .expect("Invalid structure?");

            let ltab = LoopTable::try_from((&sequence, &pairings, &emodel))
                .expect("Invalid sequence/structure combination");
            let $wname = LoopNeighbors::from((ltab, shift_policy::NoShift));
        };
    }

    fn same_moves(a: &[(Move, i32)], b: &[(Move, i32)]) -> bool {
        let sa: HashSet<_> = a.iter().copied().collect();
        let sb: HashSet<_> = b.iter().copied().collect();
        sa == sb
    }

    #[test]
    fn test_flux_after_moves() {
        let emodel = ViennaRNA::default();
        let rmodel = Metropolis::new(emodel.temperature(), 1.0, Some(1.0), Some(1.0));
        let policy = shift_policy::ThreeAndFour;

        let sequence = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").expect("Invalid sequence?");
        let pairings =     PairTable::try_from(".((((....)))).((((........))))...............").expect("Invalid structure?");

        let ltab = LoopTable::try_from((&sequence, &pairings, &emodel))
            .expect("Invalid sequence/structure combination");
        let walker = LoopNeighbors::from((ltab, policy));
 
        let mut simulator = SSA::from((walker, &rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng(), 100., |_, _, _, w| {
            steps += 1;
            println!("{}, {}", w.current_structure(), w.current_energy());
            let moves1 = w.propose_moves().collect::<Vec<_>>();
            let p = PairTable::try_from(&w.current_structure()).expect("Invalid structure?");
            let ltab = LoopTable::try_from((&sequence, &p, &emodel))
                .expect("Invalid sequence/structure combination");
            let walker = LoopNeighbors::from((ltab, policy));
            let moves2 = walker.propose_moves().collect::<Vec<_>>();
            assert!(same_moves(&moves1, &moves2));
            true
        });
        assert!(steps > 0, "Simulation must perform at least one step");
    }


    #[test]
    fn test_simple_ssa_simulation() {
        setup_ssa_input!(walker, rmodel, "GUACACGUCG", "..........");
        let mut rng = StdRng::seed_from_u64(42);
        let mut simulator = SSA::from((walker, &rmodel));

        let mut steps = 0;
        simulator.simulate(&mut rng, 1.0, |t, tinc, flux, _| {
            steps += 1;
            println!(
                "Step {}: t={:.4}, dt={:.4}, flux={:.3e}",
                steps, t, tinc, flux
            );
            true
        });
        assert!(steps > 0, "Simulation must perform at least one step");
    }

    #[test]
    fn test_cotr_ssa_simulation() {
        setup_ssa_input!(walker, rmodel, "GUACACGUCG", "......");
        let mut rng = StdRng::seed_from_u64(42);
        let mut simulator = SSA::from((walker, &rmodel));

        simulator.co_simulate(&mut rng, 
            &[4000.0, 4000.0, 4000.0, 4000.0, 4000.0], 
            |t, tinc, flux, w| {
                println!("{} {:8.2} {:14.8e} {:14.8e} {:15.8e}",
                    w,
                    w.current_energy() as f64 / 100.,
                    t,
                    tinc,
                    1.0 / flux,
                );
                true
        });
    }
}


