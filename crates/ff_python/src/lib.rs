#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::types::PyModule;
use pyo3::prelude::*;

pub mod energy_exports;
pub mod structure_exports;
pub mod kinetics_exports;

#[pymodule]
fn _fuzzyfold(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<energy_exports::Fuzzyfold>()?;
    m.add_class::<energy_exports::ViennaRNA>()?;
    m.add_function(wrap_pyfunction!(energy_exports::stems, m)?)?;
    m.add_function(wrap_pyfunction!(structure_exports::partners, m)?)?;
    m.add_function(wrap_pyfunction!(structure_exports::pairs, m)?)?;
    m.add_function(wrap_pyfunction!(structure_exports::normalize, m)?)?;
    m.add_function(wrap_pyfunction!(structure_exports::closed_regions, m)?)?;
    m.add_function(wrap_pyfunction!(structure_exports::density, m)?)?;
    m.add_class::<kinetics_exports::Simulator>()?;
    m.add_class::<kinetics_exports::Explorer>()?;
    m.add_function(wrap_pyfunction!(kinetics_exports::lin_log_times, m)?)?;
    m.add_function(wrap_pyfunction!(kinetics_exports::transcription_times, m)?)?;
    Ok(())
}


