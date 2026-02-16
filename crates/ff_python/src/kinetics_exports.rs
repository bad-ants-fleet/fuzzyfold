use ff_kinetics::RateModel;
use ff_kinetics::SSA;
use pyo3::exceptions::PyValueError;
use pyo3::types::PyModule;
use pyo3::prelude::*;

use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA as VRNA;
use ff_kinetics::Macrostate;
use ff_kinetics::Metropolis;
use ff_kinetics::LoopNeighbors;
use ff_kinetics::shift_policy::ShiftPolicy;


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
}

#[pymethods]
impl Simulator {

    #[pyo3(signature = (
        sequence,
        source,
        target,
        num_sim=100,
        t_max=1000.0
    ))]
    fn simulate(
        &self,
        py: Python<'_>,
        sequence: &str,
        source: &str,
        target: &str,
        num_sim: usize,
        t_max: f64,
    ) -> PyResult<(f64, Vec<f64>)> {

        // Convert sequence
        let seq = NucleotideVec::try_from(sequence)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // Convert source structure
        let source_pt = MultiPairTable::try_from(source)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // Convert target macrostate
        let target_macro = Macrostate::try_from(target)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        // Build move set (using stored shift_cfg)
        let moves = LoopNeighbors::try_from((
            &seq,
            &source_pt,
            &self.energy_model,
        ))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let rate_model = &self.rate_model;

        // 🔥 Release GIL during heavy computation
        let hitting_times = py.allow_threads(|| {

            let mut times = Vec::with_capacity(num_sim);

            for _ in 0..num_sim {

                let mut ssa = SSA::from((moves.clone(), rate_model));

                let t_hit = ssa.simulate_until(
                    &target_macro,
                    t_max,
                );

                times.push(t_hit);
            }

            times
        });

        let mean = if hitting_times.is_empty() {
            0.0
        } else {
            hitting_times.iter().sum::<f64>() / hitting_times.len() as f64
        };

        Ok((mean, hitting_times))
    }
}


