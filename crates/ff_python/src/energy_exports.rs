use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA as VRNA;
use ff_energy::parameters::RNA_TURNER_2004;

#[pyfunction]
pub fn eval(sequence: &str, structure: &str) -> PyResult<i32> {
    let model = VRNA::default();

    let sequence = NucleotideVec::try_from(sequence)
        .map_err(|e| PyValueError::new_err(format!("Invalid sequence: {:?}", e)))?;

    let pairings = MultiPairTable::try_from(structure)
        .map_err(|e| PyValueError::new_err(format!("Invalid structure: {:?}", e)))?;

    let energy = model.energy_of_structure(&sequence, &pairings)
        .map_err(|e| PyValueError::new_err(format!("Energy evaluation errored: {:?}", e)))?;

    Ok(energy)
}

#[pyclass]
pub struct ViennaRNA {
    model: VRNA,
}

#[pymethods]
impl ViennaRNA {

    #[new]
    #[pyo3(signature = (
        celsius=37.0,
    ))]
 
    fn new(celsius: f64) -> Self {
        Self { 
            model: VRNA::from_thermo_params(&RNA_TURNER_2004, celsius)
        }
    }

    fn energy_of_structure(
        &self,
        sequence: &str,
        structure: &str,
    ) -> PyResult<f64> {

        let sequence = NucleotideVec::try_from(sequence)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let pairings = MultiPairTable::try_from(structure)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let energy = self.model.energy_of_structure(&sequence, &pairings)
            .map_err(|e| PyValueError::new_err(format!("Energy evaluation errored: {:?}", e)))?;
 
        Ok(energy as f64 / 100.0)
    }
}


