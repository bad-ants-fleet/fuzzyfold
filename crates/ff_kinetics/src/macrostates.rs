use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io;
use std::sync::Arc;
use std::path::PathBuf;
use rustc_hash::FxHashMap;
use rand::Rng;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;

use crate::{K0, KB};

/// Represents a **macrostate**, i.e. an ensemble of secondary structures
/// sharing a common label or coarse-grained feature.
///
/// A `Macrostate` groups multiple concrete structures (`DotBracketVec`) that
/// belong to the same basin or equivalence class in the folding landscape.
/// Each entry in the ensemble is stored with:
///
/// - its **energy** value (`i32`, in dcal/mol = 100 * kcal/mol),
/// - its **probability within the macrostate** (`f64`, proportional to exp(-E/RT)).
///
/// The optional `ensemble_energy` field stores the free energy of the
/// macrostate (-kt ln(Q)) calculated at initialization. 
///
/// # Fields
/// - `name`: Identifier or label for this macrostate `α` (e.g. `"Local minimum 1"`
///   or `"Hairpin A"`).
/// - `ensemble`: Mapping from secondary structure representations `s` to `(E(s),
///   P(s|α))`.
/// - `ensemble_energy`: The free energy of the macrostate `(P(α))`.
///
/// # Notes
/// - The `MacrostateRegisty` initializes a "**catch-all** macrostate", which is
///   called "Unassigned", has no structures and ensemble_energy = None.
/// - Will panic if a structure is not well-formed.
///
/// # Example
/// ```rust
/// use ff_kinetics::Macrostate;
/// use ff_structure::DotBracketVec;
/// use ff_energy::NucleotideVec;
/// use ff_energy::ViennaRNA;
///
/// let energy_model = ViennaRNA::default();
/// let seq = NucleotideVec::try_from("CGCAAAGCG").unwrap();
/// let db1 = DotBracketVec::try_from("(((...)))").unwrap();
/// let db2 = DotBracketVec::try_from("((.....))").unwrap();
/// let db3 = DotBracketVec::try_from(".((...)).").unwrap();
///
/// let ms = Macrostate::from_list(
///    "test",
///    &seq,
///    &[db1, db2, db3],
///    &energy_model,
/// );
/// ```
#[derive(Debug)]
pub struct Macrostate {
    name: String,
    ensemble: FxHashMap<DotBracketVec, (i32, f64)>,
    ensemble_energy: Option<f64>,
}

impl Macrostate {
    /// Initialize a catch-all macrostate. 
    /// (This is the default macrostate when initializing a MacrostateRegisty.)
    pub fn new_catch_all(name: &str) -> Self {
        Macrostate { 
            name: name.to_owned(),
            ensemble: FxHashMap::default(),
            ensemble_energy: None,
        }
    }

    pub fn from_list<E: EnergyModel>(
        name: &str, 
        sequence: &NucleotideVec,
        structures: &[DotBracketVec], 
        energy_model: &E, 
    ) -> Self {
        let mut ensemble = FxHashMap::default();
        let rt = KB * (K0 + energy_model.temperature());

        let mut q_sum = 0.0;
        for dbv in structures {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = energy_model.energy_of_structure(sequence, &pt)
                .expect("Broken energy evaluation!");
            let q = (-en as f64 / 100.0 / rt).exp();
            ensemble.insert(dbv.clone(), (en, q));
            q_sum += q;
        }
        // Turn partition function contributions into probabilities.
        for (_dbv, (_en, prob)) in ensemble.iter_mut() {
            *prob /= q_sum;
        }
        Self {
            name: name.to_owned(),
            ensemble,
            ensemble_energy: Some(-rt * q_sum.ln()),
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn ensemble(&self) -> &FxHashMap<DotBracketVec, (i32, f64)> {
        &self.ensemble
    }

    pub fn ensemble_energy(&self) -> Option<f64> {
        self.ensemble_energy
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

    pub fn get_lowest_microstate(&self) -> Option<&DotBracketVec> {
        self.ensemble
            .iter()
            .min_by(|(_, (e_a, _)), (_, (e_b, _))| {
                e_a.partial_cmp(e_b).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(dbv, _)| dbv)
    }

    /// Randomly pick a structure according to its probability in the ensemble.
    pub fn get_random_microstate(&self) -> Option<DotBracketVec> {
        if self.ensemble.is_empty() {
            return None;
        }

        debug_assert!({ // Assert correct normalization
            let sum: f64 = self.ensemble.values().map(|(_, p)| p).sum();
            (sum - 1.0).abs() < 1e-6
        }, "Ensemble probabilities are not normalized!");

        let mut r = rand::rng().random::<f64>(); // random in [0, 1)

        for (dbv, &(_, p)) in &self.ensemble {
            if r < p {
                return Some(dbv.clone());
            }
            r -= p;
        }
        eprintln!("WARNING: rounding error observed. This should be rare!");
        self.ensemble.keys().next().cloned()
    }
}

/// A registy to collect macrostate definitions.
pub struct MacrostateRegistry<E: EnergyModel> {
    sequence: NucleotideVec,
    energy_model: Arc<E>,
    /// By convention: macrostates[0] = unassigned.
    macrostates: Vec<Macrostate>,
}

impl<E: EnergyModel> From<(NucleotideVec, Arc<E>)> for MacrostateRegistry<E> {
    fn from((sequence, energy_model): (NucleotideVec, Arc<E>)) -> Self {
        let macrostates = vec![Macrostate::new_catch_all("Unassigned")];

        MacrostateRegistry {
            sequence,
            energy_model,
            macrostates,
        }
    }
}

fn io_err(msg: &str, src: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{} in {}", msg, src))
}

impl<E: EnergyModel> MacrostateRegistry<E> {

    /// High-level entry: read one or more macrostate files from disk.
    pub fn insert_files(&mut self, files: &[PathBuf]) -> io::Result<()> {
        for file in files {
            let fh = File::open(file)?;
            let reader = BufReader::new(fh);
            self.insert_from_reader(reader, &file.display().to_string())?;
        }
        Ok(())
    }

    pub fn insert_from_reader<R: BufRead>(&mut self, reader: R, source: &str
    ) -> io::Result<()> {

        let mut lines = reader.lines();

        let name = {
            let line = lines
                .next()
                .ok_or_else(|| io_err("Missing first input line (header)", source))??;
            let token = line
                .trim()
                .strip_prefix('>')
                .ok_or_else(|| io_err("Header line must start with '>'", source))?
                .split_whitespace()
                .next()
                .ok_or_else(|| io_err("Header line is empty", source))?;
            if token.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
                token.to_string()
            } else {
                return Err(io_err("Header name must be alphanumeric", source));
            }
        };

        let file_seq = {
            let seq_line = lines
                .next()
                .ok_or_else(|| io_err("Missing second input line (sequence)", source))??;
            NucleotideVec::try_from(seq_line.trim())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
        };

        if file_seq != self.sequence {
            return Err(io_err("Sequence does not match input sequence", source));
        }

        let mut structures = Vec::new();
        let mut warned_trailing = false;
        for (lineno, line_raw) in lines.enumerate() {
            let line = line_raw?.trim().to_string();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let structure_str = match line.split_once(char::is_whitespace) {
                Some((s, _)) => {
                    if !warned_trailing {
                        eprintln!("Warning: trailing data after dot-bracket structures is ignored in {}.", source);
                        warned_trailing = true;
                    }
                    s
                },
                None => &line,
            };

            match DotBracketVec::try_from(structure_str) {
                Ok(dbv) => structures.push(dbv),
                Err(e) => {
                    return Err(io_err(
                        &format!("Invalid secondary structure at line {}: {:?}", lineno + 3, e),
                        source,
                    ));
                }
            }
        }

        if structures.is_empty() {
            return Err(io_err("No structures found", source));
        }

        self.macrostates.push(Macrostate::from_list(
            &name,
            &self.sequence,
            &structures,
            &*self.energy_model,
        ));
        Ok(())
    }

    /// Try to classify a structure:
    /// - Returns Some(index) if exactly one macrostate contains the structure
    /// - Returns None if no macrostate matches
    /// - Panics if more than one macrostate matches (unimplemented)
    pub fn classify(&self, structure: &DotBracketVec) -> usize {
        let mut matches = Vec::new();

        for (i, ms) in self.macrostates.iter().enumerate() {
            if ms.contains(structure) {
                matches.push(i);
            }
        }

        match matches.len() {
            0 => 0,
            1 => matches[0],
            _ => {
                // For now: raise a panic, since overlapping macrostates are ambiguous
                panic!("Structure {:?} belongs to multiple macrostates - not implemented", structure);
            }
        }
    }

    pub fn sequence(&self) -> &NucleotideVec {
        &self.sequence
    }

    pub fn macrostates(&self) -> &Vec<Macrostate> {
        &self.macrostates
    }

    /// Number of macrostates, including the catch-all unassigned macrostate.
    pub fn len(&self) -> usize {
        self.macrostates.len()
    }

    //NOTE: Useless: there is always one.
    pub fn is_empty(&self) -> bool {
        self.macrostates.is_empty()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Macrostate)> {
        self.macrostates.iter().enumerate()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::{parameters::RNA_TURNER_2004, ViennaRNA};

    #[test]
    fn test_macrostate_init() {
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

        let energy_model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
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
            "lmin=lm3_bh=3.0",
            &seq, 
            &[db1.clone(), db2, db3, db4, db5, db6, db7, db8], 
            &energy_model
        );
        println!("Macrostate '{}':", macrostate.name());
        println!("  Ensemble size: {}", macrostate.len());
        assert_eq!(macrostate.len(), 8);

        let ensemble = macrostate.ensemble().clone();
        let mut ensemble: Vec<_> = ensemble.iter().collect();
        ensemble.sort_by_key(|(_, (energy, _))| *energy);
        for (dbr, (energy, prob)) in ensemble.iter() {
            println!("  {} -> E(s) = {energy}, P(s) = {prob:.4}", dbr);
        }

        assert_eq!(macrostate.ensemble().get(&db1).unwrap().0, -390);
        assert!((macrostate.ensemble().get(&db1).unwrap().1 - 0.7669).abs() < 1e-4);
    }

    #[test]
    fn test_macrostateregistry_init_and_classify() {
        let energy_model = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();

        let mut registry = MacrostateRegistry::from((seq, Arc::new(energy_model)));

        // By convention: one unassigned macrostate
        assert_eq!(registry.len(), 1);
        assert_eq!(registry.macrostates()[0].name(), "Unassigned");

        let input = b">test
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        ";
        registry.insert_from_reader(Cursor::new(input), "manual").unwrap();
        assert_eq!(registry.len(), 2);

        // Build a test macrostate with a few structures
        let s1 = DotBracketVec::try_from(".((((....)))).((((........))))...............").unwrap();
        assert_eq!(registry.classify(&s1), 1);
        let s2 = DotBracketVec::try_from(".((((....)))).((((.(....).))))...............").unwrap();
        assert_eq!(registry.classify(&s2), 1);
        let s3 = DotBracketVec::try_from(".((((....))))..(((........)))................").unwrap();
        assert_eq!(registry.classify(&s3), 1);

        // Unknown structure: should return 0 ("Unassigned")
        let s4 = DotBracketVec::try_from("..............").unwrap();
        assert_eq!(registry.classify(&s4), 0);

        // Iteration test
        let all_names: Vec<_> = registry.iter().map(|(_, ms)| ms.name().to_string()).collect();
        assert!(all_names.contains(&"Unassigned".to_string()));
        assert!(all_names.contains(&"test".to_string()));
    }
}