#![allow(unsafe_op_in_unsafe_fn)]
use pyo3::types::PyModule;
use pyo3::prelude::*;

pub mod energy_exports;
pub mod kinetics_exports;
pub mod pseudoknot_exports;

#[pymodule]
fn fuzzyfold(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<energy_exports::ViennaRNA>()?;
    m.add_class::<kinetics_exports::Simulator>()?;
    m.add_class::<pseudoknot_exports::TreeStep>()?;
    m.add_class::<pseudoknot_exports::RegionNode>()?;
    m.add_function(wrap_pyfunction!(pseudoknot_exports::region_tree_steps, m)?)?;
    m.add_function(wrap_pyfunction!(pseudoknot_exports::region_tree, m)?)?;
    m.add_function(wrap_pyfunction!(pseudoknot_exports::pseudo_energy, m)?)?;
    m.add_function(wrap_pyfunction!(pseudoknot_exports::pseudo_energy_loops, m)?)?;
    Ok(())
}

