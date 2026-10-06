use rustc_hash::FxHashSet;
use rustc_hash::FxHashMap;
use rand::Rng;
use std::sync::Arc;
use crate::Walker;

use ff_energy::NucleotideVec;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::EnergyModel;

use crate::shift_policy;
use crate::shift_policy::ShiftPolicy;
use crate::LoopNeighbors;
use crate::Macrostate;
use crate::MacrostateRegistry;
use crate::RateModel;
use crate::enum_neighbors::ApplyMove;
use crate::macrostates::pad;

#[allow(clippy::too_many_arguments)]
fn find_neighbors<E: EnergyModel, R: RateModel>(
    dbr: &DotBracketVec,
    lss_opt: Option<&LoopNeighbors<E, shift_policy::NoShift>>,
    sequence: &NucleotideVec,
    emodel: &Arc<E>,
    rate_model: &R,
    origin: &FxHashMap<DotBracketVec, (i32, f64)>,
    visited: &mut FxHashSet<DotBracketVec>,
    neighbors: &mut FxHashMap<DotBracketVec, (i32, f64)>,
) {
    // Stop if already processed
    if !visited.insert(dbr.clone()) {
        return;
    }

    // So the outer loop does not produce 
    // all loop structures unnecessarily.
    let mut lss = match lss_opt {
        Some(existing) => existing.clone(),
        None => {
            let pairings = PairTable::try_from(dbr).unwrap();
            LoopNeighbors::try_from(
                (sequence.clone(), &pairings, emodel.clone(), shift_policy::NoShift))
                .unwrap()
        }
    };

    let mut mdbr = dbr.clone();
    for (bp_move, delta) in lss.all_moves() {
        // Let's look up whether the move is 
        // worth applying based on DotBracketVec.
        mdbr.apply_move(&bp_move);

        if origin.contains_key(&mdbr) {
            lss.apply_move(&bp_move);
            find_neighbors(&mdbr, Some(&lss), 
                sequence, emodel, rate_model, origin, visited, neighbors);
            lss.apply_move(&bp_move.inverse());
        } else {
            // Rate to step out of the macrostate => P(i|alpha) * k_{i->j}
            let rate = origin.get(dbr).unwrap().1 * rate_model.rate(&bp_move, delta);
            neighbors
                .entry(mdbr.clone())
                .and_modify(|(e, k)| {
                    debug_assert!(*e == lss.current_energy() + delta);
                    *k += rate;
                })
            .or_insert((lss.current_energy() + delta, rate));
        }
        mdbr.apply_move(&bp_move.inverse());
    }
}


#[derive(Debug)]
pub struct ExitMacrostate<'a> {
    parent_macrostate: &'a Macrostate,
    ensemble: FxHashMap<DotBracketVec, (i32, f64)>,
    k_alpha: f64,
}

impl<'a, E: EnergyModel, R: RateModel, P: ShiftPolicy> From<
(&'a Macrostate, &NucleotideVec, &Arc<E>, &R, P)> for ExitMacrostate<'a> {
    fn from((parent_macrostate, sequence, energy_model, rate_model, _shift_policy): 
        (&'a Macrostate, &NucleotideVec, &Arc<E>, &R, P)) -> Self {
        let mut visited = FxHashSet::default();
        let mut ensemble = FxHashMap::default();
        for dbr in parent_macrostate.ensemble().keys() {
            let dbr = pad(dbr, sequence.len());
            find_neighbors::<E, R>(
                &dbr,
                None,
                sequence,
                energy_model,
                rate_model,
                parent_macrostate.ensemble(),
                &mut visited,
                &mut ensemble,
            );
        }

        let mut k_alpha = 0.0;
        for (dbv, (en, k_ij)) in &ensemble {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            debug_assert_eq!(*en, energy_model.energy_of_structure(sequence, &pt).unwrap());
            k_alpha += k_ij;
        }

        ExitMacrostate {
            parent_macrostate,
            ensemble,
            k_alpha,
        }
    }
}

impl<'a> ExitMacrostate<'a> {
    pub fn parent_macrostate(&self) -> &Macrostate {
        self.parent_macrostate
    }

    pub fn ensemble(&self) -> &FxHashMap<DotBracketVec, (i32, f64)> {
        &self.ensemble
    }

    pub fn k_alpha(&self) -> f64 {
        self.k_alpha
    }
 
    /// Number of secondary structures.
    pub fn len(&self) -> usize {
        self.ensemble.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ensemble.is_empty()
    }

    /// Check if a secondary structure is contained in this macrostate.
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.ensemble.contains_key(structure)
    }

    /// Randomly pick a structure according to the exit probability.
    pub fn get_random_microstate(&self) -> Option<DotBracketVec> {
        if self.ensemble.is_empty() {
            return None
        }
                // Draw a random number in [0, total)
        let mut rng = rand::rng();
        let mut t = rng.random_range(0.0..self.k_alpha);

        // Walk through ensemble and subtract until threshold crosses 0
        for (dbv, &(_, k_ij)) in &self.ensemble {
            t -= k_ij;
            if t <= 0.0 {
                return Some(dbv.clone());
            }
        }
        eprintln!("WARNING: rounding error observed. This should be rare!");
        self.ensemble.keys().next().cloned()
    }
}


pub struct ExitMacrostateRegistry<'a, E: EnergyModel, R: RateModel> {
    parent_registry: &'a MacrostateRegistry<E>,
    rate_model: &'a R,
    /// By convention: macrostates[0] = unassigned.
    exit_macrostates: Vec<ExitMacrostate<'a>>,
}

impl<'a, E: EnergyModel, R: RateModel>
    From<(&'a MacrostateRegistry<E>, &'a R)>
    for ExitMacrostateRegistry<'a, E, R>
{
    fn from((parent_registry, rate_model): 
        (&'a MacrostateRegistry<E>, &'a R)) -> Self {
        let mut exit_macrostates = Vec::with_capacity(parent_registry.len());

        // Index 0 is unassigned, so just an empty placeholder.
        exit_macrostates.push(ExitMacrostate {
            parent_macrostate: &parent_registry.macrostates()[0].1,
            ensemble: FxHashMap::default(),
            k_alpha: 0.0,
        });
        let emodel = parent_registry.energy_model();

        // Compute neighbors for each real macrostate
        for (_, (_, ms)) in parent_registry.iter().skip(1) {
            //eprintln!("Calculating neighbors for macrostate #{i}: {}", ms.name());
            let exit_ms = ExitMacrostate::from((
                ms, 
                parent_registry.sequence(),
                &emodel.clone(),
                rate_model,
                shift_policy::NoShift,
            ));
            exit_macrostates.push(exit_ms);
        }

        ExitMacrostateRegistry {
            parent_registry, 
            rate_model,
            exit_macrostates,
        }
    }
}

impl<'a, E: EnergyModel, R: RateModel> ExitMacrostateRegistry<'a, E, R> {
    pub fn parent_registry(&self) -> &MacrostateRegistry<E> {
        self.parent_registry
    }

    pub fn rate_model(&self) -> &R {
        self.rate_model
    }

    pub fn exit_macrostates(&self) -> &Vec<ExitMacrostate<'a>> {
        &self.exit_macrostates
    }

    /// Number of exit_macrostates, including the catch-all unassigned macrostate.
    pub fn len(&self) -> usize {
        self.exit_macrostates.len()
    }

    //NOTE: Useless: there is always one.
    pub fn is_empty(&self) -> bool {
        self.exit_macrostates.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &ExitMacrostate<'_>)> {
        self.exit_macrostates.iter().enumerate()
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::ViennaRNA;
    use ff_energy::NucleotideVec;
    use crate::Arrhenius;

    #[test]
    fn test_exit_macrostate_init() {
        /*        
        >lmin=lm3_bh=3.0
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        .((((....)))).((((.(.....)))))...............
        .(((......))).((((........))))...............
        ..(((....)))..((((........))))...............
        .(((......)))..(((........)))................
        .(((.(...)))).((((........))))...............
        */        

        let energy_model = ViennaRNA::default();

        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();
        let db1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
        let db2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
        let db3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
        let db4 = DotBracketVec::try_from(".((((....)))).((((.(.....)))))...............").unwrap();
        let db5 = DotBracketVec::try_from(".(((......))).((((........))))...............").unwrap();
        let db6 = DotBracketVec::try_from("..(((....)))..((((........))))...............").unwrap();
        let db7 = DotBracketVec::try_from(".(((......)))..(((........)))................").unwrap();
        let db8 = DotBracketVec::try_from(".(((.(...)))).((((........))))...............").unwrap();

        let macrostate = Macrostate::from_list(
            "LM",
            &seq, 
            &[db1, db2, db3, db4, db5, db6, db7, db8], 
            &energy_model
        );

        let emodel = Arc::new(energy_model);

        let rate_model = Arrhenius::new(emodel.temperature(), 1.0, None, None);
        let neighbors = ExitMacrostate::from((
            &macrostate, 
            &seq, 
            &emodel, 
            &rate_model,
            shift_policy::NoShift
        ));
        println!("Neighbors '{}':", neighbors.parent_macrostate().name());
        println!("  Ensemble size: {}", neighbors.len());
        assert_eq!(neighbors.len(), 345);

        let ensemble = neighbors.ensemble().clone();
        let mut ensemble: Vec<_> = ensemble.iter().collect();
        ensemble.sort_by_key(|(_, (energy, _))| *energy);
        for (dbr, (energy, prob)) in ensemble.iter() {
            println!("  {} -> E(s) = {energy}, P(s|alpha) = {prob:.4}", dbr);
        }

    }

}


