use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyDict;

use std::fs;
use std::sync::Arc;
use std::path::PathBuf;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use ff_structure::DotBracketVec;
use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::check_sequence;
use ff_energy::Fuzzyfold;
use ff_kinetics::SSA;
use ff_kinetics::MoveTypes;
use ff_kinetics::Arrhenius;
use ff_kinetics::MoveSet;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::enum_neighbors::enumerate_neighbors;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::timecourse::Timecourse;
use ff_kinetics::ensemble::{Ensemble, State};
use crate::energy_exports::thermo_params;
use ff_kinetics::OutputTimes;
use ff_kinetics::extension_schedule;

fn resolve_energy_model(params: &str, celsius: f64) -> PyResult<(Fuzzyfold, bool)> {
    let (thermo, is_rna) = thermo_params(params)?;
    Ok((Fuzzyfold::from_thermo_params(thermo, celsius), is_rna))
}

#[pyclass]
pub struct Simulator {
    energy_model: Arc<Fuzzyfold>,
    rate_model: Arrhenius,
    is_rna: bool,
}


#[pymethods]
impl Simulator {
    #[new]
    #[pyo3(signature = (
        params = "rna_extended",
        celsius=37.0,
        k0=1e5,
        k3ws=0.0,
        k4ws=0.0,
        k_atoi=0.0,
    ))]
    fn new(
        params: &str,
        celsius: f64,
        k0: f64,
        k3ws: f64,
        k4ws: f64,
        k_atoi: f64,
    ) -> PyResult<Self> {
        let (energy_model, is_rna) = resolve_energy_model(params, celsius)?;

        if k0 < 0.0 || k3ws < 0.0 || k4ws < 0.0 || k_atoi < 0.0 {
            return Err(PyValueError::new_err(
                "Rate constants must be non-negative",
            ));
        }

        let rate_model = Arrhenius::new(celsius, k0)
            .with_three_way(k3ws)
            .with_four_way(k4ws)
            .with_a_to_i(k_atoi);

        Ok(Self {
            energy_model: Arc::new(energy_model),
            rate_model,
            is_rna,
        })
    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
    ))]
    fn simulate_trajectory(
        &self,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
    ) -> PyResult<SimulationIterator> {

        let (sequence, start_pt, num_ext) = parse_inputs(self, sequence, start)?;
        let end = t_ext.map_or(0.0, |t| t * num_ext as f64) + t_end;
        let times = extension_schedule(t_ext, num_ext, end).map_err(PyValueError::new_err)?;
        build_iterator(sequence, &start_pt, Arc::clone(&self.energy_model), self.rate_model, times)
   }
   
   #[allow(clippy::too_many_arguments)]
   #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            times=None,
            num_sims=100,
            output=None,
   ))]
   fn simulate_ensemble(
       &self,
       py: Python<'_>,
       sequence: &str,
       start: Option<&str>,
       t_ext: Option<f64>,
       times: Option<Vec<f64>>,
       num_sims: usize,
       output: Option<PathBuf>,
   ) -> PyResult<Vec<(f64, FxHashMap<String, usize>)>> {
       let (sequence, start_pt, num_ext) = parse_inputs(self, sequence, start)?;
       let output_times = output_times(times, t_ext, num_ext)?;
       let schedule = extension_schedule(t_ext, num_ext, output_times.end())
           .map_err(PyValueError::new_err)?;
       let moves = LoopNeighbors::try_from((
            sequence.clone(), &start_pt, Arc::clone(&self.energy_model), self.rate_model.move_types(),
       )).map_err(|e| PyValueError::new_err(e.to_string()))?;
       let rate_model = self.rate_model;
       let edits = rate_model.move_types().a_to_i;

       let results: Vec<Vec<State>> = py.detach(|| {
            (0..num_sims).into_par_iter().map(|_| {
                let mut ssa = SSA::from((moves.clone(), rate_model));
                ssa.sample(&mut rand::make_rng::<SmallRng>(), &schedule, &output_times,
                    |w| (w.current_sequence().clone(), w.current_structure(), w.current_energy()))
            }).collect()
       });
       let output_times = output_times.times();
       let mut ensemble = Ensemble::new(output_times);
       for trajectory in results {
           ensemble.add_trajectory(trajectory);
       }

       if let Some(output) = output {
           ensemble.write_drf(output.with_extension("drf"), edits)?;
       }

       // The structure, with edits followed by its sequence (as far as transcribed).
       let counts = ensemble.counts().iter().map(|counts| counts.iter().map(|((seq, s, _), &c)| {
           let key = if edits {
               format!("{s} {}", &seq.to_string()[..s.len()])
           } else {
               s.to_string()
           };
           (key, c)
       }).collect::<FxHashMap<String, usize>>());
       Ok(output_times.iter().copied().zip(counts).collect())
    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            times=None,
            num_sims=100,
            macrostates=vec![],
            output=None,
            num_threads=None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn simulate_timecourse(
        &self,
        py: Python<'_>,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        times: Option<Vec<f64>>,
        num_sims: usize,
        macrostates: Vec<PathBuf>,
        output: Option<PathBuf>,
        num_threads: Option<usize>,
    ) -> PyResult<Py<PyAny>> {

       let (sequence, start_pt, num_ext) = parse_inputs(self, sequence, start)?;
       let output_times = output_times(times, t_ext, num_ext)?;
       let schedule = extension_schedule(t_ext, num_ext, output_times.end())
           .map_err(PyValueError::new_err)?;
       let energy_model = Arc::clone(&self.energy_model);
       let rate_model = self.rate_model;

       let mut registry = MacrostateRegistry::from((sequence.clone(), energy_model.clone()));
       registry.insert_files(&macrostates, t_ext.is_some())
           .map_err(|e| PyValueError::new_err(e.to_string()))?;

       let pool = match num_threads {
            Some(0) => return Err(PyValueError::new_err("num_threads must be positive")),
            Some(n) => Some(rayon::ThreadPoolBuilder::new().num_threads(n).build()
                .map_err(|e| PyValueError::new_err(e.to_string()))?),
            None => None,
       };

       let moves = LoopNeighbors::try_from((
            sequence.clone(), &start_pt, Arc::clone(&energy_model), rate_model.move_types(),
       )).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let master = py.detach(|| {
            let run_one = |_: usize| -> Vec<usize> {
                let mut ssa = SSA::from((moves.clone(), rate_model));
                ssa.sample(&mut rand::make_rng::<SmallRng>(), &schedule, &output_times,
                    |w| registry.classify(w.current_sequence(), &w.current_structure()))
            };

            // Count trajectories as they finish rather than keeping them in memory.
            let empty = || Timecourse::new(output_times.times(), registry.len());
            let run_all = || (0..num_sims).into_par_iter().map(run_one)
                .fold(empty, |mut tc, trajectory| { tc.add_trajectory(&trajectory); tc })
                .reduce(empty, |mut a, b| { a.merge(&b); a });
            match &pool {
                Some(pool) => pool.install(run_all),
                None => run_all(),
            }
       });

       // Same data files as ff-timecourse: <output>.nxy, <output>.tln
       if let Some(output) = output {
            let mut nxy = Vec::new();
            master.write_nxy(&mut nxy, &registry)?;
            fs::write(output.with_extension("nxy"), nxy)?;
            master.write_tln(output.with_extension("tln"), &registry)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
       }

       // The counts of each macrostate (length, name), in registry order.
       let counts = PyDict::new(py);
       for (m, (len, ms)) in registry.iter() {
            let column: Vec<usize> = master.counts.iter().map(|row| row[m]).collect();
            counts.set_item((*len, ms.name()), column)?;
       }
       let counter: Vec<usize> = (0..master.times.len()).map(|t| master.num_trajectories(t)).collect();

       let kwargs = PyDict::new(py);
       kwargs.set_item("t_sep", output_times.t_sep())?;
       kwargs.set_item("t_ext", t_ext)?;
       kwargs.set_item("seq_length", sequence.len())?;
       kwargs.set_item("start_length", start_pt.len())?;
       let timecourse = py.import("fuzzyfold.timecourse")?.getattr("Timecourse")?
            .call((&master.times, counts, counter), Some(&kwargs))?;
       Ok(timecourse.unbind())
    }

}


#[pyclass]
pub struct Explorer {
    energy_model: Arc<Fuzzyfold>,
    is_rna: bool,
    three_way_shifts: bool,
    four_way_shifts: bool,
}

#[pymethods]
impl Explorer {
    #[new]
    #[pyo3(signature = (
        params = "rna_extended",
        celsius=37.0,
        three_way_shifts=false,
        four_way_shifts=false,
    ))]
    fn new(params: &str, celsius: f64, three_way_shifts: bool, four_way_shifts: bool) -> PyResult<Self> {
        let (energy_model, is_rna) = resolve_energy_model(params, celsius)?;
        Ok(Self {
            energy_model: Arc::new(energy_model),
            is_rna,
            three_way_shifts,
            four_way_shifts,
        })
    }

    #[pyo3(signature = (
            sequence,
            structure,
            delta=None,
            maxdist=None,
            sorted=false,
            name=None,
            output=None,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn explore(
        &self,
        py: Python<'_>,
        sequence: &str,
        structure: &str,
        delta: Option<f64>,
        maxdist: Option<usize>,
        sorted: bool,
        name: Option<String>,
        output: Option<PathBuf>,
    ) -> PyResult<Vec<(String, f64)>> {

        // Validate the macrostate name before the (possibly long) enumeration.
        let ms_name = match (&output, name) {
            (None, Some(_)) => return Err(PyValueError::new_err(
                "`name` requires `output` (the macrostate file to write)")),
            (None, None) => None,
            (Some(_), Some(n)) => Some(n),
            (Some(path), None) => Some(path.file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .ok_or_else(|| PyValueError::new_err("`output` has no file name"))?),
        };
        if let Some(n) = &ms_name
            && (n.is_empty() || !n.chars().all(|c| c.is_ascii_alphanumeric() || c == '_')) {
            return Err(PyValueError::new_err(format!(
                "Macrostate name '{}' must consist of letters, digits and '_'", n)));
        }

        let (maxdelta, maxsteps) = match (delta, maxdist) {
            (Some(d), None) => ((d * 100.0).round() as i32, usize::MAX),
            (None, Some(n)) => (i32::MAX / 2, n),
            (Some(d), Some(n)) => ((d * 100.0).round() as i32, n),
            (None, None) => {
                return Err(PyValueError::new_err(
                    "at least one of `delta` or `maxdist` must be provided",
                ));
            }
        };

        let seq = match self.is_rna {
            true => NucleotideVec::try_from_rna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            false => NucleotideVec::try_from_dna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };
        let structure_db = DotBracketVec::try_from(structure)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        check_sequence(&seq, &structure_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let pairings = MultiPairTable::try_from(&structure_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let energy_model = Arc::clone(&self.energy_model);
        let shifts = MoveTypes { three_way: self.three_way_shifts, four_way: self.four_way_shifts, a_to_i: false };

        let mut results: Vec<(String, f64)> = py.detach(move || -> Result<_, String> {
            let mut out: Vec<(String, f64)> = Vec::new();
            let mut moves = LoopNeighbors::try_from((seq, &pairings, energy_model, shifts))
                .map_err(|e| e.to_string())?;
            enumerate_neighbors(&mut moves, maxdelta, maxsteps, |db, en| {
                out.push((db.to_string(), en as f64 / 100.0));
            });
            Ok(out)
        }).map_err(PyValueError::new_err)?;

        if sorted {
            results.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        }

        if let (Some(path), Some(ms_name)) = (output, ms_name) {
            let mut content = format!(">{}\n{}\n", ms_name, sequence.trim());
            // Input structure first (it is always among the results), then the rest.
            let (root, others): (Vec<_>, Vec<_>) = results.iter()
                .partition(|(db, _)| db == structure);
            for (db, en) in root.into_iter().take(1).chain(others) {
                content.push_str(&format!("{} {:6.2}\n", db, en));
            }
            fs::write(&path, content)?;
        }
        Ok(results)
    }
}

fn parse_inputs(
        sim: &Simulator,
        sequence: &str,
        start: Option<&str>,
    ) -> PyResult<(NucleotideVec, MultiPairTable, usize)> {
        let sequence = match sim.is_rna {
            true => NucleotideVec::try_from_rna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            false => NucleotideVec::try_from_dna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        let start_db = DotBracketVec::try_from(start.unwrap_or("."))
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // The start structure covers the first part of the sequence.
        if start_db.len() > sequence.len() {
            return Err(PyValueError::new_err("start is longer than sequence"));
        }
        check_sequence(&sequence[..start_db.len()], &start_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let start_pt = MultiPairTable::try_from(&start_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let num_ext = sequence.len() - start_db.len();
        Ok((sequence, start_pt, num_ext))
}

/// The given output times, or the default ones (co-transcriptional with t_ext).
fn output_times(times: Option<Vec<f64>>, t_ext: Option<f64>, num_ext: usize) -> PyResult<OutputTimes> {
    match (times, t_ext) {
        (Some(times), _) => OutputTimes::new(times),
        (None, Some(t_ext)) => OutputTimes::transcription(t_ext, num_ext, 1.0),
        (None, None) => Ok(OutputTimes::default()),
    }.map_err(PyValueError::new_err)
}

/// Output times: 0, t_lin evenly up to t_sep, t_log in log space up to t_end.
#[pyfunction]
#[pyo3(signature = (t_sep=1e-4, t_end=1.0, t_lin=50, t_log=50))]
pub fn lin_log_times(t_sep: f64, t_end: f64, t_lin: usize, t_log: usize) -> PyResult<Vec<f64>> {
    Ok(OutputTimes::lin_log(t_sep, t_end, t_lin, t_log).map_err(PyValueError::new_err)?.times().to_vec())
}

/// Co-transcriptional output times: one per extension up to the full
/// length, then 50 in log space up to t_end after it.
#[pyfunction]
#[pyo3(signature = (t_ext, num_ext, t_end=1.0))]
pub fn transcription_times(t_ext: f64, num_ext: usize, t_end: f64) -> PyResult<Vec<f64>> {
    Ok(OutputTimes::transcription(t_ext, num_ext, t_end).map_err(PyValueError::new_err)?.times().to_vec())
}




fn build_iterator(
    seq: NucleotideVec,
    start_pt: &MultiPairTable,
    energy_model: Arc<Fuzzyfold>,
    rate_model: Arrhenius,
    times: Vec<f64>,
) -> PyResult<SimulationIterator> {
    let moves = LoopNeighbors::try_from((seq, start_pt, energy_model, rate_model.move_types()))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let ssa = SSA::from((moves, rate_model));

    Ok(SimulationIterator {
        sequence: ssa.move_set().current_sequence().clone(),
        ssa,
        rng: rand::make_rng::<SmallRng>(),
        times,
        elapsed: 0.0,
        finished: false,
    })
}

#[pyclass]
pub struct SimulationIterator {
    ssa: SSA<LoopNeighbors<Fuzzyfold>, Arrhenius>,
    rng: SmallRng,
    /// The sequence of the last state (it changes with A-to-I edits).
    sequence: NucleotideVec,
    times: Vec<f64>,
    elapsed: f64,
    finished: bool,
}

#[pymethods]
impl SimulationIterator {

    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    /// The sequence of the last state: the input sequence, with the A-to-I
    /// edits made so far (full length).
    #[getter]
    fn sequence(&self) -> String {
        self.sequence.to_string()
    }

    /// The next state: (structure, energy, time, waiting time, flux). It
    /// lasts until the next move, or until the end of the current length
    /// (the next extension, or the end of the simulation).
    fn __next__(mut slf: PyRefMut<Self>) -> Option<(String, f64, f64, f64, f64)> {
        let this: &mut Self = &mut slf;
        if this.finished {
            return None;
        }
        let rest = this.times[0];
        let (drawn, flux) = this.ssa.waiting_time(&mut this.rng).unwrap_or((f64::INFINITY, 0.0));
        let tinc = drawn.min(rest);
        let moves = this.ssa.move_set();
        if moves.current_sequence() != &this.sequence {
            this.sequence = moves.current_sequence().clone();
        }
        let state = (
            moves.to_string(),
            moves.current_energy() as f64 / 100.0, // kcal/mol
            this.elapsed,
            tinc,
            flux,
        );
        this.elapsed += tinc;
        // A move after the end of the current length is not applied.
        if drawn < rest {
            this.ssa.apply_random_move(&mut this.rng);
            this.times[0] -= tinc;
        } else {
            this.times.remove(0);
            match this.times.is_empty() {
                true => this.finished = true,
                false => this.ssa.apply_extension(),
            }
        }
        Some(state)
    }
}


