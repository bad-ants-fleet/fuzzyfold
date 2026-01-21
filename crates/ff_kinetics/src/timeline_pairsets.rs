use std::fmt;
use std::sync::Arc;
use std::error::Error;
use ff_energy::EnergyModel;
use nohash_hasher::IntMap;
use ff_structure::DotBracketVec; 

use crate::macrostates_pairsets::MacrostateRegistry;

#[derive(Debug)]
pub enum TimelineError {
    Io(std::io::Error),
    Json(serde_json::Error),
    TimepointCountMismatch { found: usize, expected: usize },
    TimeMismatch { file_time: f64, expected_time: f64 },
    MacrostateNotFound(String),
}

impl fmt::Display for TimelineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {}", e),
            Self::Json(e) => write!(f, "JSON parse error: {}", e),
            Self::TimepointCountMismatch { found, expected } =>
                write!(f, "Timeline file has {found} timepoints, expected {expected}"),
            Self::TimeMismatch { file_time, expected_time } =>
                write!(f, "Time mismatch: {file_time} vs {expected_time}"),
            Self::MacrostateNotFound(name) =>
                write!(f, "Macrostate '{name}' not found in registry"),
        }
    }
}

impl Error for TimelineError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Json(e) => Some(e),
            _ => None,
        }
    }
}


impl From<std::io::Error> for TimelineError {
    fn from(e: std::io::Error) -> Self { Self::Io(e) }
}

impl From<serde_json::Error> for TimelineError {
    fn from(e: serde_json::Error) -> Self { Self::Json(e) }
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

pub struct Timeline<'a, E: EnergyModel> {
    /// Registry of all macrostates (used to classify structures)
    pub registry: Arc<MacrostateRegistry<'a, E>>,

    /// One `Timepoint` per output time in the simulation
    pub points: Vec<Timepoint>,
}

impl<'a, E: EnergyModel> Timeline<'a, E> {
    /// Build a new empty timeline for given times and an existing macrostate registry.
    pub fn new(times: &[f64], registry: Arc<MacrostateRegistry<'a, E>>) -> Self {
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

    pub fn merge(&mut self, other: Timeline<'a, E>) {
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

impl<'a, E: EnergyModel> fmt::Display for Timeline<'a, E> {
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
                let e_a = self.registry.macrostates()[*a_idx].ensemble_energy();
                let e_b = self.registry.macrostates()[*b_idx].ensemble_energy(); 
                e_a.partial_cmp(&e_b).unwrap_or(std::cmp::Ordering::Equal)
            });

            // Sort ensemble by energy (you could make this configurable)
            for (m_idx, count) in entries {
                let occu = count as f64 / total as f64;

                let name = self.registry.macrostates()[m_idx].name();
                let energy = self.registry.macrostates()[m_idx].ensemble_energy().unwrap_or(0.0);

                writeln!(
                    f,
                    "{:13.9} {:5} {:12.8} {:>10} {:>25}",
                    time,
                    m_idx,
                    occu,
                    format!("{:10.2}", energy),
                    name,
                )?;
            }
        }
        Ok(())
    }
}

