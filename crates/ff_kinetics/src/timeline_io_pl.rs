use std::fs;
use std::result;
use std::sync::Arc;
use std::path::Path;
use serde::{Serialize, Deserialize};

use ff_energy::EnergyModel;

use crate::timeline_pl::Timepoint;
use crate::timeline_pl::Timeline;
use crate::timeline_pl::TimelineError;
use crate::macrostates_pl::MacrostateRegistryPL;

#[derive(Serialize, Deserialize)]
pub struct SerializableTimeline {
    points: Vec<SerializableTimePoint>,
}

#[derive(Serialize, Deserialize)]
pub struct SerializableTimePoint {
    time: f64,
    ensemble: Vec<(String, usize)>, // (macrostate name, count)
    counter: usize,
}

fn same_time(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

impl<E: EnergyModel> Timeline<E> {
    pub fn to_serializable(&self) -> SerializableTimeline {
        SerializableTimeline {
            points: self.points.iter().map(|tp| {
                let ensemble = tp.ensemble.iter()
                    .map(|(id, count)| {
                        let name = self.registry.macrostates()[*id].name().to_string();
                        (name, *count)
                    })
                    .collect();
                SerializableTimePoint {
                    time: tp.time,
                    ensemble,
                    counter: tp.counter,
                }
            }).collect()
        }
    }

    /// Load a timeline from a JSON file, checking against the provided registry
    pub fn from_file<P: AsRef<std::path::Path>>(
        path: P,
        times: &[f64],
        registry: Arc<MacrostateRegistryPL<E>>,
    ) -> result::Result<Self, TimelineError> {
        let data = fs::read_to_string(path)?;
        let serial: SerializableTimeline = serde_json::from_str(&data)?;

        // Sanity check: number of timepoints must match
        if serial.points.len() != times.len() {
            return Err(TimelineError::TimepointCountMismatch {
                found: serial.points.len(),
                expected: times.len(),
            });
        }

        let mut timeline = Timeline::new(times, Arc::clone(&registry));

        for (tp, serial_tp) in timeline.points.iter_mut().zip(serial.points) {
            if same_time(tp.time, serial_tp.time) {
                return Err(TimelineError::TimeMismatch {
                    file_time: serial_tp.time,
                    expected_time: tp.time,
                });
            }

            for (name, count) in serial_tp.ensemble {
                // Look up macrostate by name in registry
                if let Some((idx, _m)) = registry.iter().find(|(_, m)| m.name() == name) {
                    *tp.ensemble.entry(idx).or_insert(0) += count;
                    tp.counter += count;
                } else {
                    return Err(TimelineError::MacrostateNotFound(name));
                }
            }
        }
        Ok(timeline)
    }

    /// Load timeline data from a JSON (.tln) file, appending any timepoints not yet
    /// present in `self.points`. Macrostates are matched by name against the registry;
    /// unknown names are mapped to the Unassigned macrostate (index 0 by convention).
    ///
    /// # Notes
    /// - Call `finalize()` after all loads to ensure timepoints are sorted by time.
    /// - Can be called multiple times to accumulate data from multiple files.
    pub fn load_tln_data<P: AsRef<Path>>(
        &mut self,
        path: P,
        t_rescale: f64,
    ) -> Result<(), TimelineError> {
        let data = fs::read_to_string(path)?;
        let serial: SerializableTimeline = serde_json::from_str(&data)?;

        for serial_tp in serial.points {
            let time = serial_tp.time*t_rescale;
            let t_idx = if let Some(idx) = self.points.iter().position(|tp| same_time(tp.time, time)) {
                idx
            } else {
                self.points.push(Timepoint::new(time));
                self.points.len() - 1
            };

            for (name, count) in serial_tp.ensemble {
                let macro_idx = self.registry.iter()
                    .find(|(_, m)| m.name() == name)
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);

                let tp = &mut self.points[t_idx];
                *tp.ensemble.entry(macro_idx).or_insert(0) += count;
                tp.counter += count;
            }
        }
        Ok(())
    }

    /// Load occupancy data from an NXY-format reader, converting fractional occupancies
    /// to integer trajectory counts via `(occupancy * o_counter).round()`.
    ///
    /// # Format
    /// The first line must be a header of whitespace-separated column names, optionally
    /// prefixed with `#`. The first column must be named `time`, followed by one column
    /// per macrostate. Column names are matched against the registry by name; unknown
    /// names are mapped to the Unassigned macrostate (index 0 by convention).
    ///
    /// # Arguments
    /// - `reader`     — any `BufRead` source (file, stdin, in-memory buffer)
    /// - `t_rescale`  — multiplicative factor applied to the time column (use `1.0` for no rescaling)
    /// - `o_counter`  — total number of trajectories used to reconstruct integer counts from occupancies
    ///
    /// # Errors
    /// Returns `TimelineError` on I/O failure, float parse error, or column count mismatch.
    ///
    /// # Notes
    /// - Timepoints not yet present in `self.points` are inserted on the fly.
    /// - Call `finalize()` after all loads to ensure timepoints are sorted by time.
    /// - Can be called multiple times to accumulate data from multiple files.
    pub fn load_nxy_data<P: AsRef<Path>>(
        &mut self,
        path: P,
        t_rescale: f64,
        o_counter: usize,
    ) -> Result<(), TimelineError> {
        let content = fs::read_to_string(path)?;
        let mut map: Option<Vec<usize>> = None;
        let unassigned_idx = 0; // by convention.
        for line in content.lines() {
            if line.is_empty() {
                continue;
            }

            // Parse header
            if map.is_none() {
                let cols: Vec<&str> = line.trim_start_matches('#').split_whitespace().collect();
                // skip the "time" column
                map = Some(cols[1..].iter().map(|name| {
                    self.registry.macrostates()
                        .iter()
                        .position(|m| m.name() == *name)
                        .unwrap_or(unassigned_idx)
                }).collect());
                continue;
            }
            let map = map.as_ref().unwrap();

            let values: Vec<f64> = line
                .split_whitespace()
                .map(|s| s.parse::<f64>())
                .collect::<Result<_, _>>()?;

            let time = values[0]*t_rescale;
            let t_idx = if let Some(idx) = self.points.iter().position(|tp| same_time(tp.time, time)) {
                idx
            } else {
                self.points.push(Timepoint::new(time));
                self.points.len() - 1
            };

            for (col_idx, &occupancy) in values[1..].iter().enumerate() {
                let macro_idx = map[col_idx];
                let count = (occupancy * o_counter as f64).round() as usize;
                if count > 0 {
                    let tp = &mut self.points[t_idx];
                    *tp.ensemble.entry(macro_idx).or_insert(0) += count;
                    tp.counter += count;
                }
            }
        }
        Ok(())
    }

    pub fn finalize(&mut self) {
        self.points.sort_by(|a, b| 
            a.time.partial_cmp(&b.time).unwrap_or(std::cmp::Ordering::Equal)
        );
    }
}
