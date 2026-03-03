use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::EnergyModel;
use ff_energy::ViennaRNA as VRNA;
use ff_energy::parameters::RNA_EXTENDED;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::DNA_MATHEWS_2004;

#[pyclass]
pub struct ViennaRNA {
    model: VRNA,
    is_rna: bool,
}

#[pymethods]
impl ViennaRNA {

    #[new]
    #[pyo3(signature = (
        params = "rna_default",
        celsius=37.0,
    ))]
 
    fn new(params: &str, celsius: f64) -> PyResult<Self> {
        let mut is_rna = true;
        let thermo = match params {
            "rna_default" => &RNA_TURNER_2004,
            "rna_extended" => &RNA_EXTENDED,
            "dna" => {
                is_rna = false;
                &DNA_MATHEWS_2004
            },
            _ => {
                return Err(PyValueError::new_err(
                    format!(
                        "Unknown parameter set '{}'. \
                         Valid options are: 'rna_default', 'rna_extended', 'dna'.",
                        params
                    )
                ));
            }
        };

        Ok(Self { 
            model: VRNA::from_thermo_params(thermo, celsius),
            is_rna,
        })
    }

    fn energy_of_structure(
        &self,
        sequence: &str,
        structure: &str,
    ) -> PyResult<f64> {

        let sequence = match self.is_rna {
            true => NucleotideVec::try_from_rna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            false => NucleotideVec::try_from_dna(sequence)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        let pairings = MultiPairTable::try_from(structure)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let energy = self.model.energy_of_structure(&sequence, &pairings)
            .map_err(|e| PyValueError::new_err(format!("Energy evaluation errored: {:?}", e)))?;
 
        Ok(energy as f64 / 100.0)
    }
}


