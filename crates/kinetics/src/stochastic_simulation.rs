use std::fmt;
use rand::Rng;
use nohash_hasher::IntMap;
use energy::EnergyModel;
use crate::LoopStructure;

pub const KB: f64 = 0.001987204285; // kcal/(mol*K)
pub const K0: f64 = 273.15;

pub trait KineticModel {
    /// Given ΔE (in kcal/mol) and maybe other info, return the rate constant
    fn rate(&self, delta_e: i32) -> f64;
    fn log_rate(&self, delta_e: i32) -> f64 {
        // Default, better be overwitten.
        self.rate(delta_e).ln()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Metropolis {
    kt: f64, // k_B * T in kcal/mol
    k0: f64,
}

impl Metropolis {
    pub fn new(celsius: f64, k0: f64) -> Self {
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            kt: KB * t_kelvin,
            k0,
        }
    }
}

impl KineticModel for Metropolis {
    fn rate(&self, delta_e: i32) -> f64 {
        if delta_e <= 0 {
            self.k0
        } else {
            self.k0 * ((-delta_e as f64 / 100.) / self.kt).exp()
        }
    }

    fn log_rate(&self, delta_e: i32) -> f64 {
        if delta_e <= 0 {
            self.k0.ln()
        } else {
            self.k0.ln() + ((-delta_e as f64 / 100.) / self.kt)
        }
    }
}

fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY { return b; }
    if b == f64::NEG_INFINITY { return a; }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Compute log(exp(a) - exp(b)) safely, requires a >= b.
/// Returns -inf if the result is numerically zero or negative.
fn log_sub(a: f64, b: f64) -> Option<f64> {
    if b == f64::NEG_INFINITY {
        return Some(a);
    }
     // allow small epsilon to absorb roundoff
    if b > a + 1e-12 {
        return None; // inconsistent state, recompute needed
    }

    let gap = a - b;
    //if gap < 1e-12 {
    //    return None; // too close, cancellation risk
    //}

    let diff = (-gap).exp(); // in (0, 1]
    Some(a + (1.0 - diff).ln())
}

fn log_sum_exp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY; // empty set
    }
    m + (xs.iter().map(|&x| (x - m).exp()).sum::<f64>()).ln()
}

#[derive(Debug, Clone, PartialEq)]
pub enum Reaction {
    Add {
        i: usize,
        j: usize,
        delta_e: i32,
        log_rate: f64,
    },
    Del {
        i: usize,
        j: usize,
        delta_e: i32,
        log_rate: f64,
    },
}

impl Reaction {
    pub fn new_add<K: KineticModel>(model: &K, 
        i: usize, j: usize, delta_e: i32
) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Add { i, j, delta_e, log_rate: rate }
    }

    pub fn new_del<K: KineticModel>(model: &K, 
        i: usize, j: usize, delta_e: i32) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Del { i, j, delta_e, log_rate: rate }
    }

    pub fn ij(&self) -> (usize, usize) {
        match self {
            Reaction::Add { i, j, .. } => (*i, *j),
            Reaction::Del { i, j, .. } => (*i, *j),
        }
    }

    pub fn log_rate(&self) -> f64 {
        match self {
            Reaction::Add { log_rate, .. } => *log_rate,
            Reaction::Del { log_rate, .. } => *log_rate,
        }
    }

    pub fn delta_e(&self) -> i32 {
        match self {
            Reaction::Add { delta_e, .. } => *delta_e,
            Reaction::Del { delta_e, .. } => *delta_e,
        }
    }

}

pub struct LoopStructureSSA<'a, M: EnergyModel, K: KineticModel> {
    loopstructure: LoopStructure<'a, M>, // owns the RNA folding state
    ratemodel: &'a K,
    log_flux: f64,
    pair_flux: Option<f64>,
    loop_flux: Option<f64>,
    per_loop_flux: IntMap<usize, f64>,
    per_loop_rxns: IntMap<usize, Vec<Reaction>>,
    pair_rxns: IntMap<usize, Reaction>
}

impl<'a, M, K> fmt::Debug for LoopStructureSSA<'a, M, K>
where
    M: EnergyModel,
    K: KineticModel + fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("LoopStructureSSA")
            .field("ratemodel", &self.ratemodel)   // prints Debug for kinetic model
            .field("loopstructure", &format!("{}", self.loopstructure))
            .field("flux", &self.log_flux)
            //.field("num_reactions", &self.reactions.len())
            .finish()
    }
}

impl<'a, M: EnergyModel, K: KineticModel> From<(LoopStructure<'a, M>, &'a K)>
    for LoopStructureSSA<'a, M, K>
{
    fn from((loopstructure, ratemodel): (LoopStructure<'a, M>, &'a K)) -> Self {
        let mut per_loop_flux = IntMap::default();
        let mut per_loop_rxns = IntMap::default();
        let mut loop_logs = Vec::new();

        for (lli, add_neighbors) in loopstructure.loop_neighbors().iter() {
            let mut logs = Vec::new();
            let mut lrxns = Vec::new();
            for &(i, j, delta) in add_neighbors {
                let rxn = Reaction::new_add(ratemodel, i, j, delta);
                logs.push(rxn.log_rate());
                lrxns.push(rxn);
            }
            if lrxns.len() > 0 {
                let lflux = log_sum_exp(&logs);
                per_loop_flux.insert(*lli, lflux);
                loop_logs.push(lflux);
            }
            per_loop_rxns.insert(*lli, lrxns);
        }

        let mut pair_rxns: IntMap<usize, Reaction> = IntMap::default();
        let mut pair_logs = Vec::new();
        for (i, j, delta) in loopstructure.get_del_neighbors() {
            let rxn = Reaction::new_del(ratemodel, i, j, delta);
            pair_logs.push(rxn.log_rate());
            pair_rxns.insert(i, rxn);
        }

        let pair_flux = if pair_logs.len() > 0 {
            Some(log_sum_exp(&pair_logs))
        } else {
            None
        };

        let loop_flux = if loop_logs.len() > 0 {
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

impl<'a, M: EnergyModel, K: KineticModel> LoopStructureSSA<'a, M, K> {
    pub fn current_structure(&self) -> String {
        format!("{}", self.loopstructure)
    }

    fn recompute_flux(&mut self) {
        //println!("{}", self.current_structure());
        //println!("Recomputing flux: T:{} L:{:?} P:{:?}",
        //    self.log_flux, self.loop_flux, self.pair_flux);
        let loops: Vec<f64> = self.per_loop_flux.values().cloned().collect();
        let pairs: Vec<f64> = self.pair_rxns.values().map(|rxn| rxn.log_rate()).collect();
        self.loop_flux = if loops.len() > 0 { Some(log_sum_exp(&loops)) } else { None };
        self.pair_flux = if pairs.len() > 0 { Some(log_sum_exp(&pairs)) } else { None };
        self.log_flux = match (self.pair_flux, self.loop_flux) {
            (Some(pf), None) => pf,
            (None, Some(lf)) => lf,
            (Some(pf), Some(lf)) => log_add(pf, lf),
            _ => panic!("no flux at all?"),
        };
        //println!("Recomputed  flux: T:{} L:{:?} P:{:?}",
        //    self.log_flux, self.loop_flux, self.pair_flux);
    }
   
    pub fn remove_loop_reaction(&mut self, lli: usize) {
        let rxns = self.per_loop_rxns.remove(&lli).expect("Reaction must exist.");
        if rxns.len() == 0 {
            debug_assert!(self.per_loop_flux.remove(&lli).is_none());
            return
        }
        let lflux = self.per_loop_flux.remove(&lli)
            .expect("The lflux to be removed.");
        if self.per_loop_flux.len() > 0 {
            self.loop_flux = Some(log_sub(self.loop_flux.unwrap(), lflux).expect("lf, now that one should be fine."));
            self.log_flux = log_sub(self.log_flux, lflux).expect("tf, now that one should be fine.");
        } else {
            self.loop_flux = None;
            //NOTE: no log_flux update! Will be recomputed.
        }
    }

    pub fn remove_pair_reaction(&mut self, pli: usize) {
        let old_rxn = self.pair_rxns.remove(&pli).expect("The reaction to be removed.");
        let lrate = old_rxn.log_rate();

        if self.pair_rxns.len() > 0 {
            self.pair_flux = Some(log_sub(self.pair_flux.unwrap(), lrate).expect("pf, now that one should be fine."));
            self.log_flux = log_sub(self.log_flux, lrate).expect("tf, now that one should be fine.");
        } else {
            self.pair_flux = None;
            //NOTE: no log_flux update! Will be recomputed.
        }

        let (i, j) = old_rxn.ij();
        let &lli_outer = self.loopstructure.loop_lookup().get(&i).expect("i -> lli outer");
        let &lli_inner = self.loopstructure.loop_lookup().get(&j).expect("j -> lli inner");
        self.remove_loop_reaction(lli_inner);
        self.remove_loop_reaction(lli_outer);
    }

    pub fn insert_loop_reactions(&mut self, 
        lli: usize, 
        add_neighbors: Vec<(usize, usize, i32)>
    ) {
        let mut logs = Vec::new();
        let mut lrxns = Vec::new();
        for (i, j, delta) in add_neighbors {
            let rxn = Reaction::new_add(self.ratemodel, i, j, delta);
            logs.push(rxn.log_rate());
            lrxns.push(rxn);
        }
        if lrxns.len() > 0 {
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

    pub fn update_pair_reactions(&mut self, change: Vec<(usize, usize, i32)>) {
        for (i, j, delta) in change {
            // then it is an update, otherwise insert!
            if let Some(old) = self.pair_rxns.remove(&i) {
                let lrate = old.log_rate();
                if self.pair_rxns.len() > 0 {
                    self.pair_flux = Some(log_sub(self.pair_flux.unwrap(), lrate).expect("upf, now that one should be fine."));
                    self.log_flux = log_sub(self.log_flux, lrate).expect("utf, now that one should be fine.");
                } else {
                    self.pair_flux = None;
                    //NOTE: no log_flux update! Will be recomputed.
                }
            } 
            let rxn = Reaction::new_del(self.ratemodel, i, j, delta);
            let lrate = rxn.log_rate();
            if self.pair_rxns.len() > 0 {
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
        F: FnMut(f64, f64, f64, &LoopStructure<'a, M>),
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
            callback(t, tinc, flux, &self.loopstructure);
            t += tinc;

            // sample reaction, probably the bottleneck for now
            let log_thresh = self.log_flux + rng.random::<f64>().ln(); // ln(u) ≤ 0
            let mut acc = f64::NEG_INFINITY;
            let mut chosen = None;
           
            if let Some(pf) = self.pair_flux {
                if pf >= log_thresh {
                    for (pli, rxn) in self.pair_rxns.iter() {
                        acc = log_add(acc, rxn.log_rate());
                        if acc >= log_thresh {
                            chosen = Some((*pli, rxn.clone()));
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
                                chosen = Some((*lli, rxn.clone()));
                                break 'outer;
                            }
                        }
                    } else {
                        acc = next_acc;
                    }
                }
            }

            if let Some((idx, rxn)) = chosen {
                match rxn {
                    Reaction::Add { i, j, .. } => {
                        self.remove_loop_reaction(idx);
                        let ((lli, ami), (llj, amj), pair_changes) = self
                            .loopstructure.apply_add_move(i, j);
                        self.insert_loop_reactions(lli, ami);
                        self.insert_loop_reactions(llj, amj);
                        self.update_pair_reactions(pair_changes);
                    },
                    Reaction::Del { i, j, .. } => {
                        self.remove_pair_reaction(idx);
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

