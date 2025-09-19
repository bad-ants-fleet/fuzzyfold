use std::fmt;
use std::sync::Arc;
use ahash::AHashMap;
use nohash_hasher::IntMap;
use structure::{DotBracketVec, PairTable}; 
use energy::{EnergyModel, NucleotideVec};

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use crate::{K0, KB};

pub fn load_macrostates<M: EnergyModel>(
    files: &[PathBuf],
    sequence: &NucleotideVec,
    model: Option<&M>,
) -> MacrostateRegistry {
    let mut registry = MacrostateRegistry::new();

    for file in files {
        let fh = File::open(file)
            .unwrap_or_else(|_| panic!("Failed to open macrostate file {:?}", file));
        let mut lines = BufReader::new(fh).lines();

        // --- Parse header line (must exist)
        let header_line = lines
            .next()
            .expect("Macrostate file is missing a header line")
            .expect("Failed to read header line");

        let name = header_line
            .strip_prefix('>')
            .unwrap_or_else(|| panic!("First line must start with '>' in {:?}", file))
            .trim()
            .to_string();

        // --- Parse sequence line (must exist)
        let seq_line = lines
            .next()
            .expect("Macrostate file is missing a sequence line")
            .expect("Failed to read sequence line");

        let file_seq = NucleotideVec::from_lossy(&seq_line.trim());
        assert_eq!(
            file_seq, *sequence,
            "Sequence in macrostate file {:?} does not match provided input sequence",
            file
        );

        // --- Parse structures (must be at least one)
        let mut structures = Vec::new();
        for (lineno, line) in lines.enumerate() {
            let line = line.expect("Failed to read line");
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            match DotBracketVec::try_from(line) {
                Ok(dbv) => structures.push(dbv),
                Err(e) => panic!(
                    "Error parsing dot-bracket at line {} in {:?}: {:?}",
                    lineno + 3, // +3 to account for header+sequence lines
                    file,
                    e
                ),
            }
        }

        assert!(
            !structures.is_empty(),
            "Macrostate file {:?} does not contain any structures",
            file
        );

        // Precompute RT once if we have a model
        let mut reg = StructureRegistry::from_list(&name, &structures);

        if let Some(model) = model {
            let rt = KB * (K0 + model.temperature());
            reg.assign_energy(sequence, model, rt)
        }

        registry.insert(Macrostate::Explicit(reg));
    }

    registry
}

#[derive(Debug)]
pub enum Macrostate {
    /// Explicit registry of structures in the macrostate.
    Explicit(StructureRegistry),

    /// Placeholder for constraint-based macrostates (not yet implemented).
    Constraint(DotBracketVec),
}

impl Macrostate {
    pub fn name(&self) -> &str {
        match self {
            Macrostate::Explicit(registry) => {
                &registry.name()
            }
            Macrostate::Constraint(_) => {
                todo!("");
            }
        }
    }

    pub fn energy(&self) -> Option<f64> {
        match self {
            Macrostate::Explicit(registry) => {
                registry.energy()
            }
            Macrostate::Constraint(_) => {
                todo!("");
            }
        }
    }
}

#[derive(Debug, Default)]
pub struct StructureRegistry {
    name: String,
    pool: Vec<DotBracketVec>,
    lookup: AHashMap<DotBracketVec, usize>,
    energy: Option<f64>,
}

impl StructureRegistry {
    pub fn from_list(name: &str, structures: &[DotBracketVec]) -> Self {
        let mut pool = Vec::with_capacity(structures.len());
        let mut lookup = AHashMap::with_capacity(structures.len());

        for (i, dbv) in structures.iter().enumerate() {
            pool.push(dbv.clone());
            lookup.insert(dbv.clone(), i);
        }

        Self {
            name: name.to_owned(),
            pool,
            lookup,
            energy: None,
        }
    }

    pub fn assign_energy<M: EnergyModel>(&mut self, sequence: &NucleotideVec, model: &M, rt: f64,) {
        let mut q_sum = 0.0;
        for dbv in self.pool.iter() {
            let pt = PairTable::try_from(dbv)
                .expect("Invalid dot-bracket for energy evaluation");
            let en = model.energy_of_structure(sequence, &pt) as f64 / 100.0;
            q_sum += (-en / rt).exp();
        }
        self.energy = Some(-rt * q_sum.ln());
    }


    /// Return the name of this registry (macrostate label)
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the name of this registry (macrostate label)
    pub fn energy(&self) -> Option<f64> {
        self.energy
    }

    /// Check whether a structure is contained in this registry.
    pub fn contains(&self, structure: &DotBracketVec) -> bool {
        self.lookup.contains_key(structure)
    }

    /// Return the number of unique structures in this registry.
    pub fn len(&self) -> usize {
        self.pool.len()
    }

    /// Return an iterator over all stored structures.
    pub fn iter(&self) -> impl Iterator<Item = &DotBracketVec> {
        self.pool.iter()
    }
}

#[derive(Debug)]
pub struct MacrostateRegistry {
    macrostates: Vec<Macrostate>,
}

impl MacrostateRegistry {
    /// Start with an empty registry
    pub fn new() -> Self {
        // Start with an "unassigned" catch-all
        let unassigned = Macrostate::Explicit(StructureRegistry {
            name: "Unassigned".to_string(),
            pool: Vec::new(),
            lookup: AHashMap::new(),
            energy: None,
        });

        Self {
            macrostates: vec![unassigned],
        }
    }

    /// Insert a new macrostate
    pub fn insert(&mut self, macrostate: Macrostate) {
        self.macrostates.push(macrostate);
    }

    /// Try to classify a structure:
    /// - Returns Some(index) if exactly one macrostate contains the structure
    /// - Returns None if no macrostate matches
    /// - Panics if more than one macrostate matches (unimplemented)
    pub fn classify(&self, structure: &DotBracketVec) -> usize {
        let mut matches = Vec::new();

        for (i, ms) in self.macrostates.iter().enumerate() {
            match ms {
                Macrostate::Explicit(reg) => {
                    if reg.contains(structure) {
                        matches.push(i);
                    }
                }
                Macrostate::Constraint(_pattern) => {
                    // Not implemented yet
                }
            }
        }

        match matches.len() {
            0 => 0,
            1 => matches[0],
            _ => {
                // For now: raise a panic, since overlapping macrostates are ambiguous
                panic!("Structure {:?} belongs to multiple macrostates — not implemented", structure);
            }
        }
    }

    /// Get a macrostate by index
    pub fn get(&self, idx: usize) -> &Macrostate {
        &self.macrostates[idx]
    }

    /// Number of macrostates
    pub fn len(&self) -> usize {
        self.macrostates.len()
    }

    /// Iterate over all macrostates
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Macrostate)> {
        self.macrostates.iter().enumerate()
    }
}

/// One time point with its ensemble of macrostates.
#[derive(Debug)]
pub struct Timepoint {
    /// Absolute time in seconds
    pub time: f64,
    /// Mapping from macrostate index → number of trajectories in this state
    pub ensemble: IntMap<usize, usize>,
    /// Total number of observations recorded at this timepoint
    pub counter: usize,
}

impl Timepoint {
    /// Create a new empty timepoint
    pub fn new(time: f64) -> Self {
        Self {
            time,
            ensemble: IntMap::default(),
            counter: 0,
        }
    }

    /// Add a count for the given macrostate index
    pub fn add(&mut self, macro_idx: usize) {
        *self.ensemble.entry(macro_idx).or_insert(0) += 1;
        self.counter += 1;
    }

    /// Get the count for a specific macrostate (or 0 if not present)
    pub fn count(&self, macro_idx: usize) -> usize {
        *self.ensemble.get(&macro_idx).unwrap_or(&0)
    }

    /// Return the occupancy (fraction of total) for a macrostate
    pub fn occupancy(&self, macro_idx: usize) -> f64 {
        if self.counter == 0 {
            0.0
        } else {
            self.count(macro_idx) as f64 / self.counter as f64
        }
    }

    /// Iterate over all macrostate counts
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.ensemble.iter().map(|(k, v)| (*k, *v))
    }

}

#[derive(Debug)]
pub struct Timeline {
    /// Registry of all macrostates (used to classify structures)
    pub registry: Arc<MacrostateRegistry>,

    /// One `Timepoint` per output time in the simulation
    pub points: Vec<Timepoint>,
}

impl Timeline {
    /// Build a new empty timeline for given times and an existing macrostate registry.
    pub fn new(times: &[f64], registry: Arc<MacrostateRegistry>) -> Self {
        let points = times.iter().map(|&t| Timepoint::new(t)).collect();
        Self { registry, points }
    }

    /// Classify a structure and add it to the timeline at the given time index.
    /// Returns `Some(macro_idx)` if classified, `None` if unclassified.
    pub fn assign_structure(&mut self, t_idx: usize, structure: &DotBracketVec) {
        let m_idx = self.registry.classify(structure);
        self.points[t_idx].add(m_idx);
    }

    /// Get a reference to a timepoint by index.
    pub fn point(&self, t_idx: usize) -> &Timepoint {
        &self.points[t_idx]
    }

    /// Iterate over all timepoints with their index.
    pub fn iter(&self) -> impl Iterator<Item = (usize, &Timepoint)> {
        self.points.iter().enumerate()
    }

    pub fn merge(&mut self, other: Timeline) {
        assert!(
            Arc::ptr_eq(&self.registry, &other.registry),
            "Cannot merge timelines with different registries"
        );
        assert_eq!(self.points.len(), other.points.len(),
        "Cannot merge timelines with different numbers of timepoints");

        for (self_tp, other_tp) in self.points.iter_mut().zip(other.points) {
            for (macro_idx, count) in other_tp.iter() {
                *self_tp.ensemble.entry(macro_idx).or_insert(0) += count;
            }
            self_tp.counter += other_tp.counter;
        }
    }
}

impl fmt::Display for Timeline {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // DRF header
        writeln!(f, "{:>13} {:>5} {:>12} {:>10} {:>25}", "time", "id", "occupancy", "energy", "macrostate")?;
        for tp in self.points.iter() {
            let time = tp.time;
            let total = tp.counter.max(1);

            // Collect ensemble into a vector so we can sort it
            let mut entries: Vec<_> = tp.iter().collect();

            // Sort by energy, None last
            entries.sort_by(|(a_idx, _), (b_idx, _)| {
                let e_a = match self.registry.get(*a_idx) {
                    Macrostate::Explicit(reg) => reg.energy(),
                    Macrostate::Constraint(_) => None,
                };
                let e_b = match self.registry.get(*b_idx) {
                    Macrostate::Explicit(reg) => reg.energy(),
                    Macrostate::Constraint(_) => None,
                };

                match (e_a, e_b) {
                    (Some(a), Some(b)) => a.partial_cmp(&b).unwrap_or(std::cmp::Ordering::Equal),
                    (Some(_), None) => std::cmp::Ordering::Less,
                    (None, Some(_)) => std::cmp::Ordering::Greater,
                    (None, None) => std::cmp::Ordering::Equal,
                }
            });


            // Sort ensemble by energy (you could make this configurable)
            for (m_idx, count) in entries {
                let occu = count as f64 / total as f64;

                // Get macrostate name (or fallback)
                let (name, energy) = match self.registry.get(m_idx) {
                    Macrostate::Explicit(reg) => (reg.name(), reg.energy()),
                    Macrostate::Constraint(_) => ("Constraint", None),
                };

                writeln!(
                    f,
                    "{:13.9} {:5} {:12.8} {:>10} {:>25}",
                    time,
                    m_idx,
                    occu,
                    energy.map_or("N/A".to_string(), |e| format!("{:10.2}", e)),
                    name,
                )?;
            }
        }
        Ok(())
    }
}


