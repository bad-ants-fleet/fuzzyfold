use ff_structure::DotBracketVec;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use rand::rng;

use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA as VRNA;
use ff_kinetics::SSA;
use ff_kinetics::shift_policy;
use ff_kinetics::Metropolis;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;


#[pyclass]
pub struct Simulator {
    energy_model: VRNA,
    rate_model: Metropolis,
}

#[pymethods]
impl Simulator {

    #[new]
    #[pyo3(signature = (
        temperature=37.0,
        k0=1e5,
        k3ws=0.0,
        k4ws=0.0,
    ))]
    fn new(
        temperature: f64,
        k0: f64,
        k3ws: f64,
        k4ws: f64,
    ) -> PyResult<Self> {

        if k0 < 0.0 || k3ws < 0.0 || k4ws < 0.0 {
            return Err(PyValueError::new_err(
                "Rate constants must be non-negative",
            ));
        }

        let mut energy_model = VRNA::default();
        energy_model.set_temperature(temperature);

        let rate_model = Metropolis::new(
            temperature,
            k0,
            Some(k3ws),
            Some(k4ws),
        );

        Ok(Self {
            energy_model,
            rate_model,
        })
    }

    #[pyo3(signature = (
        sequence,
        start,
        t_ext=None,
        t_end=1.0,
        silent=false,
    ))]
    fn simulate(
        &self,
        py: Python<'_>,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
        silent: bool,
    ) -> PyResult<Vec<(String, i32, f64, f64, f64)>> {

        // Convert sequence
        let seq = NucleotideVec::try_from(sequence)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let start = match start {
            Some(s) => DotBracketVec::try_from(s)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => DotBracketVec::try_from(".")
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        if start.len() < seq.len() && t_ext.is_none() {
            return Err(PyValueError::new_err(
                    "t_ext must be provided when the start stucture is shorter than the sequence length!",
            ));
        }

        let times: Vec<f64> = if let Some(t_ext) = t_ext {
            let mut v = vec![t_ext; sequence.len() - start.len()];
            v.push(t_end);
            v
        } else { 
            vec![t_end] 
        };

        // Convert source structure
        let start = MultiPairTable::try_from(&start)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // Build move set (using stored shift_cfg)
        let moves = LoopNeighbors::try_from((
            &seq,
            &start,
            &self.energy_model,
            shift_policy::NoShift,
        ))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let mut results = Vec::new();

        py.allow_threads(|| {
            let mut ssa = SSA::from((moves, &self.rate_model));

            if silent {
                ssa.co_simulate(&mut rng(), &times, |_, _, _, _| true);
                let cs = ssa.current_structure().to_string();
                let en = ssa.current_energy();
                results.push((cs, en, t_end, 0.0, 0.0));
            } else {
                ssa.co_simulate(&mut rng(), &times, 
                    |t, tinc, flux, w| {
                        results.push((
                                w.to_string(),
                                w.current_energy(),
                                t,
                                tinc,
                                1.0 / flux,
                        ));
                        true 
                    });
            }
        });

        Ok(results)
    }
}


