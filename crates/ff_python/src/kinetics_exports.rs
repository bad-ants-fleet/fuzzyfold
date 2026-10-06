use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use pyo3::types::{PyDict, PyList};

use std::fs;
use std::ffi::CString;
use std::sync::Arc;
use std::path::PathBuf;
use rand::SeedableRng;
use rand::rngs::SmallRng;
use rayon::prelude::*;
use rustc_hash::FxHashMap;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA;
use ff_kinetics::SSA;
use ff_kinetics::shift_policy;
use ff_kinetics::Arrhenius;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::MacrostateRegistry;
use ff_kinetics::timeline::Timeline;
use ff_energy::parameters::RNA_EXTENDED;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::DNA_MATHEWS_2004;
use ff_shared::kinetics_parsers::TimelineParameters;

fn resolve_energy_model(params: &str, celsius: f64) -> PyResult<(ViennaRNA, bool)> {
    let mut is_rna = true;
    let thermo = match params {
        "rna_turner_2004" => &RNA_TURNER_2004,
        "rna_extended" => &RNA_EXTENDED,
        "dna" => {
            is_rna = false;
            &DNA_MATHEWS_2004
        },
        _ => {
            return Err(PyValueError::new_err(
                format!(
                    "Unknown parameter set '{}'. \
                     Valid options are: 'rna_turner_2004', 'rna_extended', 'dna'.",
                    params
                )
            ));
        }
    };
    Ok((ViennaRNA::from_thermo_params(thermo, celsius), is_rna))
}

#[pyclass]
pub struct Simulator {
    energy_model: Arc<ViennaRNA>,
    rate_model: Arrhenius,
    is_rna: bool,
}


#[pymethods]
impl Simulator {
    #[new]
    #[pyo3(signature = (
        params = "rna_turner_2004",
        celsius=37.0,
        k0=1e5,
        k3ws=0.0,
        k4ws=0.0,
    ))]
    fn new(
        params: &str,
        celsius: f64,
        k0: f64,
        k3ws: f64,
        k4ws: f64,
    ) -> PyResult<Self> {
        let (energy_model, is_rna) = resolve_energy_model(params, celsius)?;

        if k0 < 0.0 || k3ws < 0.0 || k4ws < 0.0 {
            return Err(PyValueError::new_err(
                "Rate constants must be non-negative",
            ));
        }

        let rate_model = Arrhenius::new(
            celsius,
            k0,
            Some(k3ws),
            Some(k4ws),
        );

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

        let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;
        
        match (self.rate_model.k3ws().is_some(), self.rate_model.k4ws().is_some()) {
            (false, false) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::NoShift,
                SSAKind::NoShift,
            ),

            (true, false) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeWayOnly,
                SSAKind::ThreeWayOnly,
            ),

            (false, true) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::FourWayOnly,
                SSAKind::FourWayOnly,
            ),

            (true, true) => build_iterator(
                sequence,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeAndFour,
                SSAKind::ThreeAndFour,
            ),
        }
   }
   
   #[allow(clippy::too_many_arguments)]
   #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
            t_lin=None,
            t_log=50,
            t_sep=None,
            num_sims=100,
            output=None,
   ))]
   fn simulate_ensemble(
       &self,
       py: Python<'_>,
       sequence: &str,
       start: Option<&str>,
       t_ext: Option<f64>,
       t_end: f64,
       t_lin: Option<usize>,
       t_log: usize,
       t_sep: Option<f64>,
       num_sims: usize,
       output: Option<PathBuf>,
   ) -> PyResult<Vec<(f64, FxHashMap<String, usize>)>> {
       let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;

       let k3ws = self.rate_model.k3ws().is_some();
       let k4ws = self.rate_model.k4ws().is_some();
       let energy_model = Arc::clone(&self.energy_model);
       let rate_model = self.rate_model;

       let mut tl_params = TimelineParameters {
           t_ext,
           t_end,
           t_sep,
           t_lin,
           t_log,
       };

       let num_ext = sequence.len() - start_pt.len();
       let k0 = self.rate_model.k0().ok_or_else(|| PyValueError::new_err("rate model has no k0 set"))?;

       tl_params.validate(k0, num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let output_times = tl_params.get_output_times(num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let results: Result<Vec<Vec<(String, i32)>>, String> = py.detach(|| {
            let run_one = |_: usize| -> Result<Vec<(String, i32)>, String> {
                let mut structures: Vec<(String, i32)> = Vec::new();

                macro_rules! run_with_policy {
                    ($policy:expr) => {{
                        let walker = LoopNeighbors::try_from((
                            sequence.clone(), &start_pt, Arc::clone(&energy_model), $policy,
                        )).map_err(|e| format!("{:?}", e))?;

                        let mut ssa = SSA::from((walker, rate_model));
                        let mut rng = SmallRng::from_os_rng();
                        let mut t_idx = 0; 

                        ssa.co_simulate(&mut rng, &times, |t, tinc, _flux, w| {
                            while t_idx < output_times.len() && t + tinc >= output_times[t_idx] {
                                structures.push((w.to_string(), w.current_energy()));
                                t_idx += 1;
                            }
                            true
                        });
                    }};

                }

                match (k3ws, k4ws) {
                    (false, false) => run_with_policy!(shift_policy::NoShift),
                    (true,  false) => run_with_policy!(shift_policy::ThreeWayOnly),
                    (false, true)  => run_with_policy!(shift_policy::FourWayOnly),
                    (true,  true)  => run_with_policy!(shift_policy::ThreeAndFour),
                }
                Ok(structures)
            };

            (0..num_sims).into_par_iter().map(run_one).collect()
       });
       
       let results = results.map_err(PyValueError::new_err)?;

       let mut counts: Vec<FxHashMap<String, usize>> = (0..output_times.len()).map(|_| FxHashMap::default()).collect();
       // Energy per structure (energy in dcal/mol), for the .drf file.
       let mut energies: FxHashMap<String, i32> = FxHashMap::default();

       for structures in results {
           for (t_idx, (structure, energy)) in structures.into_iter().enumerate() {
               *counts[t_idx].entry(structure.clone()).or_insert(0) += 1;
               energies.entry(structure).or_insert(energy);
           }
       }

       if let Some(output) = output {
           write_drf(&output.with_extension("drf"), &output_times, &counts, &energies,
           num_sims, sequence.len())?;
       }

       Ok(output_times.iter().copied().zip(counts).collect())
    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
            t_lin=None,
            t_log=50,
            t_sep=None,
            num_sims=100,
            macrostates=vec![],
            output=None,
            num_threads=None,
    ))]
    fn simulate_timecourse(
        &self,
        py: Python<'_>,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
        t_lin: Option<usize>,
        t_log: usize,
        t_sep: Option<f64>,
        num_sims: usize,
        macrostates: Vec<PathBuf>,
        output: Option<PathBuf>,
        num_threads: Option<usize>,
    ) -> PyResult<Vec<(f64, FxHashMap<String, f64>)>> {

       let (sequence, start_pt, times) = parse_inputs(self, sequence, start, t_ext, t_end)?;

       let k3ws = self.rate_model.k3ws().is_some();
       let k4ws = self.rate_model.k4ws().is_some();
       let energy_model = Arc::clone(&self.energy_model);
       let rate_model = self.rate_model;

       let mut tl_params = TimelineParameters {
           t_ext,
           t_end,
           t_sep,
           t_lin,
           t_log,
       };

       let num_ext = sequence.len() - start_pt.len();
       let k0 = self.rate_model.k0().ok_or_else(|| PyValueError::new_err("rate model has no k0 set"))?;

       tl_params.validate(k0, num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let output_times = tl_params.get_output_times(num_ext).map_err(|e| PyValueError::new_err(e.to_string()))?;

       let mut ms = MacrostateRegistry::from((sequence.clone(), energy_model.clone()));
       insert_macrostate_files(&mut ms, &macrostates, t_ext.is_some())?;
       let registry = Arc::new(ms);

       let pool = match num_threads {
            Some(0) => return Err(PyValueError::new_err("num_threads must be positive")),
            Some(n) => Some(rayon::ThreadPoolBuilder::new().num_threads(n).build()
                .map_err(|e| PyValueError::new_err(e.to_string()))?),
            None => None,
       };

       let merged: Option<Result<Timeline<ViennaRNA>, String>> = py.detach(|| {
            let run_one = |_: usize| -> Result<Timeline<ViennaRNA>, String> {
                let thread_registry = Arc::clone(&registry);
                let mut timeline = Timeline::new(&output_times, thread_registry);

                macro_rules! run_with_policy {
                    ($policy:expr) => {{
                        let walker = LoopNeighbors::try_from((
                            sequence.clone(), &start_pt, Arc::clone(&energy_model), $policy,
                        )).map_err(|e| format!("{:?}", e))?;

                        let mut ssa = SSA::from((walker, rate_model));
                        let mut rng = SmallRng::from_os_rng();
                        let mut t_idx = 0; 

                        ssa.co_simulate(&mut rng, &times, |t, tinc, _flux, w| {
                            while t_idx < output_times.len() && t + tinc >= output_times[t_idx] {
                                let structure = w.current_structure();
                                timeline.assign_structure(t_idx, &structure);
                                t_idx += 1;
                            }
                            true
                        });
                    }};

                }

                match (k3ws, k4ws) {
                    (false, false) => run_with_policy!(shift_policy::NoShift),
                    (true,  false) => run_with_policy!(shift_policy::ThreeWayOnly),
                    (false, true)  => run_with_policy!(shift_policy::FourWayOnly),
                    (true,  true)  => run_with_policy!(shift_policy::ThreeAndFour),
                }
                Ok(timeline)
            };

            // Merge timelines as they finish rather than keeping all of them in memory.
            let run_all = || (0..num_sims).into_par_iter().map(run_one)
                .try_reduce_with(|mut a, b| { a.merge(b); Ok(a) });
            match &pool {
                Some(pool) => pool.install(run_all),
                None => run_all(),
            }
       });

       let master = match merged {
            Some(timeline) => timeline.map_err(PyValueError::new_err)?,
            None => Timeline::new(&output_times, Arc::clone(&registry)),
       };

       // Same data files as ff-timecourse: <output>.nxy, <output>.tln
       if let Some(output) = output {
            write_timecourse_files(py, &master, &output)?;
       }

       Ok(timeline_to_occupancy(&master))
    }


    #[staticmethod]
    #[pyo3(signature = (inputs, output=None))]
    fn aggregate_timecourse_results(
        py: Python<'_>,
        inputs: Vec<PathBuf>,
        output: Option<PathBuf>,
    ) -> PyResult<Vec<(f64, FxHashMap<String, f64>)>> {
        let merged = merge_tln_files(py, &inputs)?;
        let timeline = counts_to_timeline(tln_to_counts(&merged))?;

        if let Some(output) = output {
            write_tln(py, &merged, &output.with_extension("tln"))?;
            fs::write(output.with_extension("nxy"), format!("{}", timeline))?;
        }
        Ok(timeline_to_occupancy(&timeline))
    }
}


#[pyclass]
pub struct Explorer {
    energy_model: Arc<ViennaRNA>,
    is_rna: bool,
    three_way_shifts: bool,
    four_way_shifts: bool,
}

#[pymethods]
impl Explorer {
    #[new]
    #[pyo3(signature = (
        params = "rna_turner_2004",
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
        if let Some(n) = &ms_name {
            if n.is_empty() || !n.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
                return Err(PyValueError::new_err(format!(
                    "Macrostate name '{}' must consist of letters, digits and '_'", n)));
            }
        }

        let (maxdelta, maxsteps) = match (delta, maxdist) {
            (Some(d), None) => ((d * 100.0) as i32, usize::MAX),
            (None, Some(n)) => (i32::MAX / 2, n),
            (Some(d), Some(n)) => ((d * 100.0) as i32, n),
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
        let pairings = PairTable::try_from(&structure_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let energy_model = Arc::clone(&self.energy_model);
        let three_way_shifts = self.three_way_shifts;
        let four_way_shifts = self.four_way_shifts;

        let mut results: Vec<(String, f64)> = py.detach(move || -> Result<_, String> {
            let mut out: Vec<(String, f64)> = Vec::new();

            macro_rules! run_with_policy {
                ($policy:expr) => {{
                    let mut moves = LoopNeighbors::try_from(
                        (seq, &pairings, energy_model, $policy)
                    ).map_err(|e| format!("{:?}", e))?;
                    moves.generate_neighbors(maxdelta, maxsteps, |db, en| {
                        out.push((db.to_string(), en as f64 / 100.0));
                    });
                }};
            }

            match (three_way_shifts, four_way_shifts) {
                (false, false) => run_with_policy!(shift_policy::NoShift),
                (true, false)  => run_with_policy!(shift_policy::ThreeWayOnly),
                (false, true)  => run_with_policy!(shift_policy::FourWayOnly),
                (true, true)   => run_with_policy!(shift_policy::ThreeAndFour),
            }

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
        t_ext: Option<f64>,
        t_end: f64,
    ) -> PyResult<(NucleotideVec, PairTable, Vec<f64>)> {
        let sequence = match sim.is_rna {
            true => NucleotideVec::try_from_rna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            false => NucleotideVec::try_from_dna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        let start_db = match start {
            Some(s) => DotBracketVec::try_from(s)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => DotBracketVec::try_from(".")
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        if start_db.len() < sequence.len() && t_ext.is_none() {
            return Err(PyValueError::new_err(
                    "t_ext must be provided when start is shorter than sequence",
            ));
        }

        let times = if let Some(dt) = t_ext {
            let mut v = vec![dt; sequence.len() - start_db.len()];
            v.push(t_end);
            v
        } else {
            vec![t_end]
        };

        let start_pt = PairTable::try_from(&start_db)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        Ok((sequence, start_pt, times))
}


fn insert_macrostate_files(
    registry: &mut MacrostateRegistry<ViennaRNA>,
    files: &[PathBuf],
    cotrans: bool,
) -> PyResult<()> {
    for file in files {
        let text = fs::read_to_string(file)
            .map_err(|e| PyValueError::new_err(format!("{}: {}", file.display(), e)))?;
        // Keep header and sequence lines as they are, only strip structure lines.
        let cleaned: String = text.lines().enumerate()
            .map(|(i, line)| match i {
                0 | 1 => line,
                _ => line.split_whitespace().next().unwrap_or(""),
            })
            .collect::<Vec<_>>()
            .join("\n");
        let source = file.display().to_string();
        let mut num_to_remove = 0;
        loop {
            let inserted = registry.insert_from_reader(cleaned.as_bytes(), &source, num_to_remove)
                .map_err(|e| PyValueError::new_err(e.to_string()))?;
            if !inserted || !cotrans {
                break;
            }
            num_to_remove += 1;
        }
    }
    Ok(())
}


fn write_drf(
    path: &PathBuf,
    times: &[f64],
    counts: &[FxHashMap<String, usize>],
    energies: &FxHashMap<String, i32>,
    num_sims: usize,
    seqlen: usize,
) -> PyResult<()> {
    use std::io::Write;

    let padded = |s: &str| format!("{}{}", s, ".".repeat(seqlen.saturating_sub(s.len())));

    // Per time point: (structure, energy, count), sorted by energy (then structure).
    let sorted: Vec<Vec<(&String, i32, usize)>> = counts.iter().map(|ensemble| {
        let mut entries: Vec<(&String, i32, usize)> = ensemble.iter()
            .map(|(s, &c)| (s, energies[s], c))
            .collect();
        entries.sort_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(b.0)));
        entries
    }).collect();

    let mut idict: FxHashMap<String, usize> = FxHashMap::default();
    for entries in &sorted {
        for (s, _, _) in entries {
            let next_id = idict.len();
            idict.entry(padded(s)).or_insert(next_id);
        }
    }

    let mut writer = std::io::BufWriter::new(fs::File::create(path)?);
    writeln!(writer, "id time occupancy structure energy")?;
    for (t, entries) in times.iter().zip(&sorted) {
        for (s, en, count) in entries {
            writeln!(
                writer,
                "{:5} {:03.3} {:5} {} {:6.2}",
                idict[&padded(s)], t, *count as f64 / num_sims as f64, s, *en as f64 / 100.0
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn write_timecourse_files(
    py: Python<'_>,
    timeline: &Timeline<ViennaRNA>,
    output: &PathBuf,
) -> PyResult<()> {
    let tln_path = output.with_extension("tln");
    let nxy_path = output.with_extension("nxy");

    fs::write(&nxy_path, format!("{}", timeline))?;

    let macrostates = timeline.registry.macrostates();
    let points: Vec<TlnPoint> = timeline.points.iter().map(|tp| TlnPoint {
        time: tp.time,
        counter: tp.counter,
        ensemble: tp.ensemble.iter()
            .map(|(id, count)| ((macrostates[*id].0, macrostates[*id].1.name().to_string()), *count))
            .collect(),
    }).collect();
    write_tln(py, &points, &tln_path)

}

#[derive(FromPyObject)]
enum TimecourseInput {
    Path(PathBuf),
    Occupancy(Vec<(f64, FxHashMap<String, f64>)>),
}

type TimecourseCounts = (Vec<String>, Vec<(f64, usize, FxHashMap<String, usize>)>);

fn ordered_names<'a>(names: impl Iterator<Item = &'a String>) -> Vec<String> {
    let mut names: Vec<String> = names.cloned().collect();
    names.sort();
    names.dedup();
    names.sort_by_key(|n| n != "Unassigned");
    names
}

struct TlnPoint {
    time: f64,
    counter: usize,
    ensemble: FxHashMap<(usize, String), usize>,
}

fn read_tln_points(py: Python<'_>, path: &PathBuf) -> PyResult<Vec<TlnPoint>> {
    let err = |msg: &str| PyValueError::new_err(format!("{}: {}", path.display(), msg));
    let text = fs::read_to_string(path)
        .map_err(|e| err(&e.to_string()))?;
    let serial = py.import("json")?.call_method1("loads", (text,))?;

    let mut points = Vec::new();
    for tp in serial.get_item("points")?.extract::<Vec<Bound<'_, PyAny>>>()? {
        let mut ensemble: FxHashMap<(usize, String), usize> = FxHashMap::default();
        for entry in tp.get_item("ensemble")?.extract::<Vec<Vec<Bound<'_, PyAny>>>>()? {
            if entry.len() != 3 {
                return Err(err("ensemble entries must be [length, name, count]"));
            }
            *ensemble.entry((entry[0].extract()?, entry[1].extract()?)).or_insert(0)
                += entry[2].extract::<usize>()?;
        }
        points.push(TlnPoint {
            time: tp.get_item("time")?.extract()?,
            counter: tp.get_item("counter")?.extract()?,
            ensemble,
        });
    }
    Ok(points)
}

fn write_tln(py: Python<'_>, points: &[TlnPoint], path: &PathBuf) -> PyResult<()> {
    let serial_points = PyList::empty(py);
    for tp in points {
        let ensemble: Vec<(usize, &str, usize)> = tp.ensemble.iter()
            .map(|((len, name), count)| (*len, name.as_str(), *count))
            .collect();
        let point = PyDict::new(py);
        point.set_item("time", tp.time)?;
        point.set_item("ensemble", ensemble)?;
        point.set_item("counter", tp.counter)?;
        serial_points.append(point)?;
    }
    let serial = PyDict::new(py);
    serial.set_item("points", serial_points)?;
    let kwargs = PyDict::new(py);
    kwargs.set_item("indent", 2)?;
    let json: String = py.import("json")?
        .call_method("dumps", (serial,), Some(&kwargs))?
        .extract()?;
    fs::write(path, json)?;
    Ok(())
}

fn tln_to_counts(points: &[TlnPoint]) -> TimecourseCounts {
    let points: Vec<(f64, usize, FxHashMap<String, usize>)> = points.iter().map(|tp| {
        let mut counts: FxHashMap<String, usize> = FxHashMap::default();
        for ((_, name), count) in &tp.ensemble {
            *counts.entry(name.clone()).or_insert(0) += count;
        }
        (tp.time, tp.counter, counts)
    }).collect();
    let names = ordered_names(points.iter().flat_map(|(_, _, c)| c.keys()));
    (names, points)
}

fn merge_tln_files(py: Python<'_>, inputs: &[PathBuf]) -> PyResult<Vec<TlnPoint>> {
    let (first, rest) = inputs.split_first()
        .ok_or_else(|| PyValueError::new_err("No input files given"))?;
    let mut merged = read_tln_points(py, first)?;
    for path in rest {
        let points = read_tln_points(py, path)?;
        if points.len() != merged.len() {
            return Err(PyValueError::new_err(format!(
                "{} has {} time points, {} has {}",
                path.display(), points.len(), first.display(), merged.len())));
        }
        for (m, tp) in merged.iter_mut().zip(points) {
            if (m.time - tp.time).abs() > 1e-9 * m.time.abs().max(1.0) {
                return Err(PyValueError::new_err(format!(
                    "{}: time {} does not match time {} in {} (different simulation parameters?)",
                    path.display(), tp.time, m.time, first.display())));
            }
            m.counter += tp.counter;
            for (key, count) in tp.ensemble {
                *m.ensemble.entry(key).or_insert(0) += count;
            }
        }
    }
    Ok(merged)
}

fn counts_to_timeline((names, points): TimecourseCounts) -> PyResult<Timeline<ViennaRNA>> {
    let (energy_model, _) = resolve_energy_model("rna_default", 37.0)?;
    let dummy_seq = NucleotideVec::try_from("A")
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mut registry = MacrostateRegistry::from((dummy_seq, Arc::new(energy_model)));
    for name in names.iter().filter(|n| *n != "Unassigned") {
        let definition = format!(">{}\nA\n.\n", name);
        registry.insert_from_reader(definition.as_bytes(), name, 0)
            .map_err(|e| PyValueError::new_err(format!("Invalid macrostate name '{}': {}", name, e)))?;
    }
    let index: FxHashMap<String, usize> = registry.macrostates().iter().enumerate()
        .map(|(i, (_, ms))| (ms.name().to_string(), i))
        .collect();

    let times: Vec<f64> = points.iter().map(|(t, _, _)| *t).collect();
    let mut timeline = Timeline::new(&times, Arc::new(registry));
    for (tp, (_, counter, counts)) in timeline.points.iter_mut().zip(points) {
        tp.counter = counter;
        for (name, count) in counts {
            if count > 0 {
                tp.ensemble.insert(index[&name], count);
            }
        }
    }
    Ok(timeline)
}

fn infer_t_sep(times: &[f64]) -> Option<f64> {
    if times.len() < 3 {
        return None;
    }
    let step = times[1] - times[0];
    let mut i = 1;
    while i + 1 < times.len() && ((times[i + 1] - times[i]) - step).abs() <= 1e-3 * step {
        i += 1;
    }
    (i + 1 < times.len()).then(|| times[i])
}

const PLOT_ANXY_SRC: &str = include_str!("../../../examples/py-utils/plot_anxy.py");

fn occupancy_to_anxy(occupancy: &[(f64, FxHashMap<String, f64>)]) -> String {
    let names = ordered_names(occupancy.iter().flat_map(|(_, o)| o.keys()));
    let mut text = format!("time {}\n", names.join(" "));
    for (t, occ) in occupancy {
        text.push_str(&format!("{:e}", t));
        for name in &names {
            text.push_str(&format!(" {:e}", occ.get(name).copied().unwrap_or(0.0)));
        }
        text.push('\n');
    }
    text
}

fn counts_to_occupancy((_, points): TimecourseCounts) -> Vec<(f64, FxHashMap<String, f64>)> {
    points.into_iter().map(|(t, counter, counts)| {
        let occ = counts.into_iter()
            .map(|(name, c)| (name, if counter > 0 { c as f64 / counter as f64 } else { 0.0 }))
            .collect();
        (t, occ)
    }).collect()
}

fn anxy_times(text: &str) -> PyResult<Vec<f64>> {
    text.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .skip(1) // header
        .map(|l| l.split_whitespace().next().unwrap_or("").parse::<f64>()
            .map_err(|e| PyValueError::new_err(format!("invalid time value '{}': {}", l, e))))
        .collect()
}

#[pyclass]
pub struct Plotter;

#[pymethods]
impl Plotter {
    #[new]
    fn new() -> Self {
        Plotter
    }
    #[staticmethod]
    #[pyo3(signature = (
            data,
            output=None,
            title=None,
            t_sep=None,
            show=true,
            **kwargs,
    ))]
    fn plot_timecourse(
        py: Python<'_>,
        data: TimecourseInput,
        output: Option<PathBuf>,
        title: Option<String>,
        t_sep: Option<f64>,
        show: bool,
        kwargs: Option<&Bound<'_, PyDict>>,
    ) -> PyResult<Option<String>> {
        let (basename, mut formats, written) = match &output {
            None if !show => return Err(PyValueError::new_err(
                "Nothing to do: give an `output` file and/or set show=True")),
            None => (String::new(), vec![], None),
            Some(path) => {
                let ext = path.extension().map_or("pdf".to_string(), |e| e.to_string_lossy().to_string());
                let base = path.with_extension("");
                let written = format!("{}.{}", base.display(), ext);
                (base.to_string_lossy().to_string(), vec![ext], Some(written))
            },
        };
        if show {
            formats.push("show".to_string());
        }

        let anxy = match data {
            TimecourseInput::Occupancy(occupancy) => occupancy_to_anxy(&occupancy),
            TimecourseInput::Path(path) => match path.extension().and_then(|e| e.to_str()) {
                Some("tln") => occupancy_to_anxy(&counts_to_occupancy(
                    tln_to_counts(&read_tln_points(py, &path)?))),
                Some("nxy") => fs::read_to_string(&path)
                    .map_err(|e| PyValueError::new_err(format!("{}: {}", path.display(), e)))?,
                _ => return Err(PyValueError::new_err(format!(
                    "{}: expected a .tln or .nxy file", path.display()))),
            },
        };

        let t_sep = match t_sep {
            Some(t) => t,
            None => infer_t_sep(&anxy_times(&anxy)?).ok_or_else(|| PyValueError::new_err(
                "Could not infer t_sep (all time points are evenly spaced?), please provide it"))?,
        };

        let code = CString::new(PLOT_ANXY_SRC).expect("plot_anxy.py contains no NUL bytes");
        let module = PyModule::from_code(py, &code, c"plot_anxy.py", c"fuzzyfold_plot_anxy")
            .map_err(|e| {
                if e.is_instance_of::<pyo3::exceptions::PyImportError>(py) {
                    pyo3::exceptions::PyImportError::new_err(format!(
                        "Plotter needs matplotlib and numpy (pip install matplotlib): {}", e))
                } else { e }
            })?;

        let call_kwargs = PyDict::new(py);
        call_kwargs.set_item("basename", basename)?;
        call_kwargs.set_item("formats", formats)?;
        call_kwargs.set_item("title", title.unwrap_or_default())?;
        call_kwargs.set_item("t_split", t_sep)?;
        if let Some(extra) = kwargs {
            for (key, value) in extra.iter() {
                let key: String = key.extract()?;
                if matches!(key.as_str(), "stream" | "basename" | "formats" | "t_split") {
                    return Err(PyValueError::new_err(format!(
                        "'{}' is set by plot_timecourse (use data/output/t_sep instead)", key)));
                }
                call_kwargs.set_item(key, value)?;
            }
        }

        let stream = py.import("io")?.getattr("StringIO")?.call1((anxy,))?;
        module.getattr("plot_anxy")?.call((stream,), Some(&call_kwargs))?;
        Ok(written)
    }
}

fn timeline_to_occupancy(timeline: &Timeline<ViennaRNA>) -> Vec<(f64, FxHashMap<String, f64>)> {
    
    let macrostates = timeline.registry.macrostates();
    let mut name_order: Vec<&str> = Vec::new();
    let mut name_to_indices: FxHashMap<&str, Vec<usize>> = FxHashMap::default();
    for (idx, (_len, ms)) in macrostates.iter().enumerate() {
        let name = ms.name();
        name_to_indices
            .entry(name)
            .or_insert_with(|| {
                name_order.push(name);
                Vec::new()
            })
        .push(idx);
    }
    
    timeline.points.iter().map(|tp| {
        let mut occupancy: FxHashMap<String, f64> = FxHashMap::default();
        for name in &name_order {
            let indices = &name_to_indices[name];
            let count: usize = indices.iter()
                .map(|&i| tp.ensemble.get(&i).copied().unwrap_or(0))
                .sum();
            let occu = if tp.counter > 0 {
                count as f64 / tp.counter as f64
            } else { 0.0 };
            occupancy.insert(name.to_string(), occu);
        }
        (tp.time, occupancy)
    }).collect()
}

fn build_iterator<P>(
    seq: NucleotideVec,
    start_pt: &PairTable,
    energy_model: Arc<ViennaRNA>,
    rate_model: Arrhenius,
    times: Vec<f64>,
    policy: P,
    wrap: fn(SSA<LoopNeighbors<ViennaRNA, P>, Arrhenius>) -> SSAKind,
) -> PyResult<SimulationIterator>
where
    P: shift_policy::ShiftPolicy,
{
    let walker = LoopNeighbors::try_from((
        seq,
        start_pt,
        energy_model,
        policy,
    ))
    .map_err(|e| PyValueError::new_err(e.to_string()))?;

    let ssa = wrap(SSA::from((walker, rate_model)));

    Ok(SimulationIterator {
        ssa,
        rng: SmallRng::from_os_rng(),
        times,
        elapsed: 0.0,
        finished: false,
    })
}

enum SSAKind {
    NoShift(SSA<LoopNeighbors<ViennaRNA, shift_policy::NoShift>, Arrhenius>),
    ThreeWayOnly(SSA<LoopNeighbors<ViennaRNA, shift_policy::ThreeWayOnly>, Arrhenius>),
    FourWayOnly(SSA<LoopNeighbors<ViennaRNA, shift_policy::FourWayOnly>, Arrhenius>),
    ThreeAndFour(SSA<LoopNeighbors<ViennaRNA, shift_policy::ThreeAndFour>, Arrhenius>),
}

#[pyclass]
pub struct SimulationIterator {
    ssa: SSAKind,
    rng: SmallRng,
    times: Vec<f64>,
    elapsed: f64,
    finished: bool,
}

#[pymethods]
impl SimulationIterator {

    fn __iter__(slf: PyRef<Self>) -> PyRef<Self> {
        slf
    }

    fn __next__(
        mut slf: PyRefMut<Self>
    ) -> Option<(String, i32, f64, f64, f64)> {

        let this: &mut Self = &mut slf;

        if this.finished {
            return None;
        }

        let mut produced: Option<(String, i32, f64, f64, f64)> = None;

        let rng = &mut this.rng;
        let mut mytinc = 0.0;
        let mut first_pass = true;

        macro_rules! dispatch_ssa {
            ($ssa:expr) => {{
                $ssa.co_simulate(
                    rng,
                    &this.times,
                    |t, tinc, flux, w| {
                        if first_pass {
                            mytinc = tinc.min(this.times[0]);

                            produced = Some((
                                    w.to_string(),
                                    w.current_energy(),
                                    this.elapsed + t,
                                    mytinc,
                                    flux,
                            ));

                            this.elapsed += mytinc;
                            first_pass = false;
                            true
                        } else {
                            false
                        }
                    },
                    );
            }};
        }

        match &mut this.ssa {
            SSAKind::NoShift(ssa) => dispatch_ssa!(ssa),
            SSAKind::ThreeWayOnly(ssa) => dispatch_ssa!(ssa),
            SSAKind::FourWayOnly(ssa) => dispatch_ssa!(ssa),
            SSAKind::ThreeAndFour(ssa) => dispatch_ssa!(ssa),
        }

        if (this.times[0] - mytinc).abs() < f64::EPSILON {
            this.times.remove(0); 
            if this.times.is_empty() {
                this.finished = true;
            }
        } else {
            assert!(this.times[0] > mytinc);
            this.times[0] -= mytinc;
        }
        produced
    }
}


