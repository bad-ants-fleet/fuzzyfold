use pyo3::exceptions::PyValueError;
use pyo3::types::PyModule;
use pyo3::prelude::*;

use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA as VRNA;

#[pyfunction]
fn eval(sequence: &str, structure: &str) -> PyResult<i32> {
    let model = VRNA::default();

    let sequence = NucleotideVec::try_from(sequence)
        .map_err(|e| PyValueError::new_err(format!("Invalid sequence: {:?}", e)))?;

    let pairings = MultiPairTable::try_from(structure)
        .map_err(|e| PyValueError::new_err(format!("Invalid structure: {:?}", e)))?;

    Ok(model.energy_of_structure(&sequence, &pairings))
}

#[pyclass]
struct ViennaRNA {
    model: VRNA,
}

#[pymethods]
impl ViennaRNA {

    #[new]
    fn new(temp_celsius: Option<f64>) -> Self {
        let mut model = VRNA::default();
        if let Some(t) = temp_celsius {
            model.set_temperature(t);
        }
        Self { model }
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

        let energy = self.model.energy_of_structure(&sequence, &pairings);

        Ok(energy as f64 / 100.0)
    }
}

#[pymodule]
fn fuzzyfold(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<ViennaRNA>()?;
    m.add_function(wrap_pyfunction!(eval, m)?)?;
    Ok(())
}


/*
#[cfg(test)]
mod tests {
    use super::*;

    fn it_works() {
    }
}
*/

