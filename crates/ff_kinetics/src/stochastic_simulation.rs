use std::fmt;
use rand::Rng;
use nohash_hasher::IntMap;
use ff_structure::NAIDX;
use ff_energy::EnergyModel;

use crate::reaction::Reaction;
use crate::LoopStructure;
use crate::RateModel;

fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Compute log(exp(a) - exp(b)) safely, requires a >= b.
/// Returns None if the result is numerically zero or negative.
fn log_sub(a: f64, b: f64) -> Option<f64> {
    if b == f64::NEG_INFINITY {
        return Some(a);
    }
    if b > a + 1e-12 {
        return None; 
    }
    let diff = (b - a).exp(); // in (0, 1]
    Some(a + (1.0 - diff).ln())
}

fn log_sum_exp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY; // empty set
    }
    m + (xs.iter().map(|&x| (x - m).exp()).sum::<f64>()).ln()
}

pub struct LoopStructureSSA<'a, M: EnergyModel, K: RateModel> {
    loopstructure: LoopStructure<'a, M>, // owns the RNA folding state
    ratemodel: &'a K,
    log_flux: f64,
    pair_flux: Option<f64>,
    loop_flux: Option<f64>,
    per_loop_flux: IntMap<usize, f64>,
    per_loop_rxns: IntMap<usize, Vec<Reaction>>,
    pair_rxns: IntMap<NAIDX, Reaction>
}

impl<'a, M, K> fmt::Debug for LoopStructureSSA<'a, M, K>
where
    M: EnergyModel,
    K: RateModel + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructureSSA")
            .field("ratemodel", &self.ratemodel)   // prints Debug for RateModel
            .field("loopstructure", &format!("{}", self.loopstructure))
            .field("flux", &self.log_flux)
            //.field("num_reactions", &self.reactions.len())
            .finish()
    }
}

impl<'a, M: EnergyModel, K: RateModel> From<(LoopStructure<'a, M>, &'a K)>
    for LoopStructureSSA<'a, M, K>
{
    fn from((loopstructure, ratemodel): (LoopStructure<'a, M>, &'a K)) -> Self {
        let mut per_loop_flux = IntMap::default();
        let mut per_loop_rxns = IntMap::default();
        let mut loop_logs = Vec::new();

        for (lli, add_neighbors) in loopstructure.get_add_neighbors_per_loop().iter() {
            let mut logs = Vec::new();
            let mut lrxns = Vec::new();
            for &(i, j, delta) in add_neighbors {
                let rxn = Reaction::new_add(ratemodel, i, j, delta);
                logs.push(rxn.log_rate());
                lrxns.push(rxn);
            }
            if !lrxns.is_empty() {
                let lflux = log_sum_exp(&logs);
                per_loop_flux.insert(*lli, lflux);
                loop_logs.push(lflux);
            }
            per_loop_rxns.insert(*lli, lrxns);
        }

        let mut pair_rxns: IntMap<NAIDX, Reaction> = IntMap::default();
        let mut pair_logs = Vec::new();
        for (i, j, delta) in loopstructure.get_del_neighbors() {
            let rxn = Reaction::new_del(ratemodel, i, j, delta);
            pair_logs.push(rxn.log_rate());
            pair_rxns.insert(i, rxn);
        }

        let pair_flux = if !pair_logs.is_empty() {
            Some(log_sum_exp(&pair_logs))
        } else {
            None
        };

        let loop_flux = if !loop_logs.is_empty() {
            Some(log_sum_exp(&loop_logs))
        } else {
            None
        };

        let log_flux = match (pair_flux, loop_flux) {
            (Some(pf), None) => pf,
            (None, Some(lf)) => lf,
            (Some(pf), Some(lf)) => log_add(pf, lf),
            _ => panic!("no flux at all?"),
        };

        Self {
            ratemodel,
            loopstructure,
            log_flux,
            pair_flux,
            loop_flux,
            per_loop_flux,
            per_loop_rxns,
            pair_rxns,
        }
    }
}

impl<'a, M: EnergyModel, K: RateModel> LoopStructureSSA<'a, M, K> {
    pub fn current_structure(&self) -> String {
        format!("{}", self.loopstructure)
    }

    fn recompute_flux(&mut self) {
        //println!("{}", self.current_structure());
        //println!("Recomputing flux: T:{} L:{:?} P:{:?}",
        //    self.log_flux, self.loop_flux, self.pair_flux);
        let loops: Vec<f64> = self.per_loop_flux.values().cloned().collect();
        let pairs: Vec<f64> = self.pair_rxns.values().map(|rxn| rxn.log_rate()).collect();
        self.loop_flux = if !loops.is_empty() { Some(log_sum_exp(&loops)) } else { None };
        self.pair_flux = if !pairs.is_empty() { Some(log_sum_exp(&pairs)) } else { None };
        self.log_flux = match (self.pair_flux, self.loop_flux) {
            (Some(pf), None) => pf,
            (None, Some(lf)) => lf,
            (Some(pf), Some(lf)) => log_add(pf, lf),
            _ => panic!("no flux at all?"),
        };
        //println!("Recomputed  flux: T:{} L:{:?} P:{:?}",
        //    self.log_flux, self.loop_flux, self.pair_flux);
    }
   
    pub fn remove_loop_reaction(&mut self, i: NAIDX) {
        let lli = self.loopstructure.loop_lookup().get(&i).unwrap();
        let rxns = self.per_loop_rxns.remove(lli).expect("Reaction must exist.");
        if rxns.is_empty() {
            debug_assert!(self.per_loop_flux.remove(lli).is_none());
            return
        }
        let lflux = self.per_loop_flux.remove(lli)
            .expect("The lflux to be removed.");
        if !self.per_loop_flux.is_empty() {
            self.loop_flux = Some(log_sub(self.loop_flux.unwrap(), lflux).expect("lf, now that one should be fine."));
            self.log_flux = log_sub(self.log_flux, lflux).expect("tf, now that one should be fine.");
        } else {
            self.loop_flux = None;
            //NOTE: no log_flux update! Will be recomputed.
        }
    }

    pub fn remove_pair_reaction(&mut self, i: NAIDX) {
        let old_rxn = self.pair_rxns.remove(&i).expect("The reaction to be removed.");
        let lrate = old_rxn.log_rate();

        if !self.pair_rxns.is_empty() {
            self.pair_flux = Some(log_sub(self.pair_flux.unwrap(), lrate).expect("pf, now that one should be fine."));
            self.log_flux = log_sub(self.log_flux, lrate).expect("tf, now that one should be fine.");
        } else {
            self.pair_flux = None;
            //NOTE: no log_flux update! Will be recomputed.
        }
        let (i, j) = old_rxn.ij();
        self.remove_loop_reaction(i);
        self.remove_loop_reaction(j);
    }

    pub fn insert_loop_reactions(&mut self, 
        lli: usize, 
        add_neighbors: Vec<(NAIDX, NAIDX, i32)>
    ) {
        let mut logs = Vec::with_capacity(add_neighbors.len());
        let mut lrxns = Vec::with_capacity(add_neighbors.len());
        for (i, j, delta) in add_neighbors {
            let rxn = Reaction::new_add(self.ratemodel, i, j, delta);
            logs.push(rxn.log_rate());
            lrxns.push(rxn);
        }
        if !lrxns.is_empty() {
            let lflux = log_sum_exp(&logs);
            self.per_loop_flux.insert(lli, lflux);
            if self.loop_flux.is_some() {
                self.loop_flux = Some(log_add(self.loop_flux.unwrap(), lflux));
            } else {
                self.loop_flux = Some(lflux);
            }
            self.log_flux = log_add(self.log_flux, lflux);
        }
        self.per_loop_rxns.insert(lli, lrxns);
    }

    pub fn update_pair_reactions(&mut self, change: Vec<(NAIDX, NAIDX, i32)>) {
        for (i, j, delta) in change {
            // then it is an update, otherwise insert!
            if let Some(old) = self.pair_rxns.remove(&i) {
                let lrate = old.log_rate();
                if !self.pair_rxns.is_empty() {
                    self.pair_flux = Some(log_sub(self.pair_flux.unwrap(), lrate).expect("upf, now that one should be fine."));
                    self.log_flux = log_sub(self.log_flux, lrate).expect("utf, now that one should be fine.");
                } else {
                    self.pair_flux = None;
                    //NOTE: no log_flux update! Will be recomputed.
                }
            } 
            let rxn = Reaction::new_del(self.ratemodel, i, j, delta);
            let lrate = rxn.log_rate();
            if !self.pair_rxns.is_empty() {
                self.pair_flux = Some(log_add(self.pair_flux.unwrap(), lrate));
            } else {
                self.pair_flux = Some(lrate);
            }
            self.log_flux = log_add(self.log_flux, lrate);
            self.pair_rxns.insert(i, rxn);
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
        F: FnMut(f64, f64, f64, &LoopStructure<'a, M>) -> bool,
    {
        let mut t = 0.;

        while t < t_max {
            if let (Some(pf), Some(lf)) = (self.pair_flux, self.loop_flux) {
                if (log_add(pf, lf) - self.log_flux).abs() > 1e-8 {
                    self.recompute_flux();
                }
            } else { self.recompute_flux(); };

            let flux = self.log_flux.exp();
            // sample waiting time ~ Exp(flux)
            let tinc = -rng.random::<f64>().ln() / flux;

            // Callback bewore applying the waiting time.
            // If callback return's false, then abort the simulation!
            if !callback(t, tinc, flux, &self.loopstructure) {
                break;
            }

            t += tinc;

            // sample reaction, probably the bottleneck for now
            let log_thresh = self.log_flux + rng.random::<f64>().ln(); // ln(u) ≤ 0
            let mut acc = f64::NEG_INFINITY;
            let mut chosen = None;
           
            if let Some(pf) = self.pair_flux {
                if pf >= log_thresh {
                    for rxn in self.pair_rxns.values() {
                        acc = log_add(acc, rxn.log_rate());
                        if acc >= log_thresh {
                            chosen = Some(rxn.clone());
                            break;
                        }
                    }
                } else {
                    acc = pf;
                }
            }
            if chosen.is_none() {
                'outer: for (lli, lflux) in self.per_loop_flux.iter() {
                    let rxns = &self.per_loop_rxns[lli];
                    let next_acc = log_add(acc, *lflux);
                    if next_acc > log_thresh {
                        for rxn in rxns {
                            acc = log_add(acc, rxn.log_rate());
                            if acc >= log_thresh {
                                chosen = Some(rxn.clone());
                                break 'outer;
                            }
                        }
                    } else {
                        acc = next_acc;
                    }
                }
            }

            if let Some(rxn) = chosen {
                match rxn {
                    Reaction::Add { i, j, .. } => {
                        self.remove_loop_reaction(rxn.ij().0);
                        let ((lli, ami), (llj, amj), pair_changes) = self
                            .loopstructure.apply_add_move(i, j);
                        self.insert_loop_reactions(lli, ami);
                        self.insert_loop_reactions(llj, amj);
                        self.update_pair_reactions(pair_changes);
                    },
                    Reaction::Del { i, j, .. } => {
                        self.remove_pair_reaction(rxn.ij().0);
                        let ((lli, neighbors), pair_changes) = self
                            .loopstructure.apply_del_move(i, j);
                        self.insert_loop_reactions(lli, neighbors);
                        self.update_pair_reactions(pair_changes);
                    },
                }
            } else {
                panic!("No reaction chosen despite positive flux");
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
        assert!(simulator.log_flux.is_finite(), "Flux must remain finite");
    }
}



