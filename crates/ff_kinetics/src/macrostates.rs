use std::io::BufRead;
use std::io;
use std::sync::Arc;
use std::path::PathBuf;
use rustc_hash::FxHashMap;
use rand::Rng;

use ff_structure::DotBracketVec;
use ff_structure::DotBracket;
use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;

use crate::{ZERO_CELSIUS, KB};

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
/// - `name`: Identifier or label for this macrostate `A` (e.g. `"Local minimum 1"`
///   or `"Hairpin A"`).
/// - `ensemble`: Mapping from secondary structure representations `s` to `(E(s),
///   P(s|A))`.
/// - `ensemble_energy`: The free energy of the macrostate `(P(A))`.
///
/// # Notes
/// - The `MacrostateRegistry` initializes a "**catch-all** macrostate", which is
///   called "Unassigned", has no structures and ensemble_energy = None.
/// - Will panic if a structure is not well-formed.
///
/// # Example
/// ```rust
/// use ff_kinetics::Macrostate;
/// use ff_structure::DotBracketVec;
/// use ff_energy::NucleotideVec;
/// use ff_energy::Fuzzyfold;
///
/// let energy_model = Fuzzyfold::default();
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

fn truncate(dbv: DotBracketVec, num: usize) -> Option<DotBracketVec> {
    let mut truncated = dbv.to_vec();
    let mut removed = 0;
    while removed < num && truncated.last() == Some(&DotBracket::UNPAIRED) {
        truncated.pop();
        removed += 1;
    }
    if removed == num {
        return Some(DotBracketVec(truncated))
    };
    None
}

pub fn pad(dbv: &DotBracketVec, len: usize) -> DotBracketVec {
    let mut padded = dbv.to_vec();
    let extra = len.saturating_sub(padded.len());
    padded.extend(std::iter::repeat_n(DotBracket::UNPAIRED, extra));
    DotBracketVec(padded)
}


impl Macrostate {
    /// Initialize a catch-all macrostate. 
    /// (This is the default macrostate when initializing a MacrostateRegistry.)
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
        let rt = KB * (ZERO_CELSIUS + energy_model.temperature());

        let mut q_sum = 0.0;
        for dbv in structures {
            let pt = MultiPairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = energy_model.energy_of_structure(sequence, &pt)
                .expect("Broken energy evaluation!");
            let q = (-en as f64 / 100.0 / rt).exp();
            ensemble.insert(dbv.clone(), (en, q));
            q_sum += q;
        }
        // Turn partition function contributions into probabilities.
        for (_en, prob) in ensemble.values_mut() {
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
        // May speed up co-transcriptional scenarios.
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

/// A registry to collect macrostate definitions.
pub struct MacrostateRegistry<E: EnergyModel> {
    sequence: NucleotideVec,
    energy_model: Arc<E>,
    /// By convention: macrostates[0] = unassigned.
    macrostates: Vec<(usize, Macrostate)>,
}

impl<E: EnergyModel> From<(NucleotideVec, Arc<E>)> for MacrostateRegistry<E> {
    fn from((sequence, energy_model): (NucleotideVec, Arc<E>)) -> Self {
        let macrostates = vec![(sequence.len(), Macrostate::new_catch_all("Unassigned"))];

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

    /// Read macrostate files. Co-transcriptional: also insert each
    /// macrostate at every shorter length (its structures cut at the 3'
    /// end), as long as any structure remains.
    pub fn insert_files(&mut self, files: &[PathBuf], cotrans: bool) -> io::Result<()> {
        for file in files {
            let source = file.display().to_string();
            let text = std::fs::read_to_string(file).map_err(|e| io_err(&e.to_string(), &source))?;
            let mut remove = 0;
            while self.insert_from_reader(text.as_bytes(), &source, remove)? && cotrans {
                remove += 1;
            }
        }
        Ok(())
    }

    pub fn insert_from_reader<R: BufRead>(&mut self, 
        reader: R, 
        source: &str, 
        remove: usize,
    ) -> io::Result<bool> {

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
        for (lineno, line_raw) in lines.enumerate() {
            let line = line_raw?.trim().to_string();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            // Anything after the structure (e.g. its energy) is ignored.
            let structure_str = line.split_whitespace().next().unwrap_or("");

            match DotBracketVec::try_from(structure_str) {
                Ok(dbv) => { 
                    if let Some(dbv) = truncate(dbv, remove) {
                        structures.push(dbv)
                    };
                }
                Err(e) => {
                    return Err(io_err(
                        &format!("Invalid secondary structure at line {}: {:?}", lineno + 3, e),
                        source,
                    ));
                }
            }
        }

        if structures.is_empty() && remove == 0 {
            return Err(io_err("No structures found", source));
        } else if structures.is_empty() {
            return Ok(false)
        }

        // Macrostates of the same length are disjoint.
        let len = file_seq.len() - remove;
        for (l, ms) in &self.macrostates {
            if let Some(s) = structures.iter().find(|s| *l == len && ms.contains(s)) {
                return Err(io_err(&format!("Structure {} is also in macrostate {}", s, ms.name()), source));
            }
        }
        self.macrostates.push((len, Macrostate::from_list(
            &name,
            &self.sequence,
            &structures,
            &*self.energy_model,
        )));
        Ok(true)
    }

    /// The index of the macrostate that contains a structure, 0 (unassigned)
    /// if none does. Macrostates of the same length are disjoint.
    pub fn classify(&self, structure: &DotBracketVec) -> usize {
        self.macrostates.iter()
            .position(|(len, ms)| structure.len() == *len && ms.contains(structure))
            .unwrap_or(0)
    }

    pub fn sequence(&self) -> &NucleotideVec {
        &self.sequence
    }

    pub fn energy_model(&self) -> &Arc<E> {
        &self.energy_model
    }


    pub fn macrostates(&self) -> &Vec<(usize, Macrostate)> {
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
    pub fn iter(&self) -> impl Iterator<Item = (usize, &(usize, Macrostate))> {
        self.macrostates.iter().enumerate()
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;
    use ff_energy::{parameters::RNA_TURNER_2004, Fuzzyfold};

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

        let energy_model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
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
        assert_eq!(macrostate.len(), 8);
        let total: f64 = macrostate.ensemble().values().map(|&(_, p)| p).sum();
        assert!((total - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_macrostateregistry_init_and_classify() {
        let energy_model = Fuzzyfold::from_thermo_params(&RNA_TURNER_2004, 37.0);
        let seq = NucleotideVec::try_from("UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC").unwrap();

        let mut registry = MacrostateRegistry::from((seq, Arc::new(energy_model)));

        // By convention: one unassigned macrostate
        assert_eq!(registry.len(), 1);
        assert_eq!(registry.macrostates()[0].1.name(), "Unassigned");

        let input = b">test
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....)))).((((........))))...............
        .((((....)))).((((.(....).))))...............
        .((((....))))..(((........)))................
        ";
        registry.insert_from_reader(Cursor::new(input), "manual", 0).unwrap();
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

        // A structure in two macrostates is an error.
        let overlap = b">other
        UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
        .((((....))))..(((........)))................
        ";
        let err = registry.insert_from_reader(Cursor::new(overlap), "manual", 0).unwrap_err();
        assert!(err.to_string().contains("also in macrostate test"), "{err}");
        assert_eq!(registry.len(), 2);

        // Iteration test
        let all_names: Vec<_> = registry.iter().map(|(_, (_, ms))| ms.name().to_string()).collect();
        assert!(all_names.contains(&"Unassigned".to_string()));
        assert!(all_names.contains(&"test".to_string()));
    }
}
