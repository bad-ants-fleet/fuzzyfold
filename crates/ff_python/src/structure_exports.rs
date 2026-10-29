//! Secondary structures as dot-bracket strings: pairs, normal form, and the
//! closed regions of pseudoknots. Indices are positions in the string, so a
//! strand break ('+' or '&') takes one position.

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use ff_structure::ClosedRegions;
use ff_structure::DotBracketVec;
use ff_structure::MultiPairTable;
use ff_structure::MultiStruct;

fn error(e: ff_structure::StructureError) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn table(structure: &str) -> PyResult<MultiPairTable> {
    MultiPairTable::try_from(structure).map_err(error)
}

/// The partner of each position: its index, or None if it is unpaired or a
/// strand break.
#[pyfunction]
pub fn partners(structure: &str) -> PyResult<Vec<Option<usize>>> {
    Ok(table(structure)?.iter()
        .map(|s| match s {
            MultiStruct::Paired(j) => Some(j),
            _ => None,
        })
        .collect())
}

/// The pairs (i, j) with i < j, sorted by i.
#[pyfunction]
pub fn pairs(structure: &str) -> PyResult<Vec<(usize, usize)>> {
    Ok(partners(structure)?.into_iter().enumerate()
        .filter_map(|(i, p)| p.filter(|&j| i < j).map(|j| (i, j)))
        .collect())
}

/// The structure with its brackets assigned again, from '(' upward, and
/// strand breaks written as '+'. Equal structures give equal strings.
#[pyfunction]
pub fn normalize(structure: &str) -> PyResult<String> {
    Ok(DotBracketVec::try_from(&table(structure)?).map_err(error)?.to_string())
}

/// The closed regions as (i, j, nested), sorted
/// by i. A region is nested if i pairs with j, pseudoknotted otherwise.
#[pyfunction]
pub fn closed_regions(structure: &str) -> PyResult<Vec<(usize, usize, bool)>> {
    Ok(ClosedRegions::from(&table(structure)?).iter()
        .map(|r| (r.span.0, r.span.1, r.nested))
        .collect())
}

/// The pseudoknot density: the largest number of bands of one
/// pseudoknotted region that span a common position. 1 for nested
/// structures, 2 for H-types and kissing hairpins. Bands are split at strand
/// breaks, so it can depend on the order of the strands.
#[pyfunction]
pub fn density(structure: &str) -> PyResult<usize> {
    let mpt = table(structure)?;
    Ok(ClosedRegions::from(&mpt).density())
}
