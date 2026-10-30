use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use ff_structure::DotBracketVec;
use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;
use ff_energy::check_sequence;
use ff_energy::EnergyError;
use ff_energy::EnergyModel;
use ff_energy::Gap;
use ff_energy::Ring;
use ff_energy::RingDecomposition;
use ff_energy::parameters::RNAThermoParams;
use ff_energy::parameters::RNA_EXTENDED;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::DNA_MATHEWS_2004;

/// The energy model of fuzzyfold: ViennaRNA loops and pseudoknots.
#[pyclass]
pub struct Fuzzyfold {
    model: ff_energy::Fuzzyfold,
    is_rna: bool,
}

#[pymethods]
impl Fuzzyfold {

    #[new]
    #[pyo3(signature = (
        params = "rna_extended",
        celsius = 37.0,
    ))]
    fn new(params: &str, celsius: f64) -> PyResult<Self> {
        let (thermo, is_rna) = thermo_params(params)?;
        if celsius != 37.0 {
            Python::attach(|py| PyErr::warn(py, &py.get_type::<pyo3::exceptions::PyUserWarning>(),
                c"Pseudoknot parameters are given at 37 C and are not rescaled with temperature.", 1))?;
        }
        Ok(Self { model: ff_energy::Fuzzyfold::from_thermo_params(thermo, celsius), is_rna })
    }

    fn energy_of_structure(&self, sequence: &str, structure: &str) -> PyResult<f64> {
        let (sequence, pairings) = inputs(self.is_rna, sequence, structure)?;
        kcal(self.model.energy_of_structure(&sequence, &pairings))
    }

    /// The energy of each loop and pseudoknot (level 0), and for a 
    /// pseudoknot also its terms (levels 1 and 2), as 
    /// `(level, description, energy in kcal/mol)`. Level 0 sums to the 
    /// energy of the structure. Same content as `ff-eval -v`.
    fn energy_terms(&self, sequence: &str, structure: &str) -> PyResult<Vec<(usize, String, f64)>> {
        let (sequence, pairings) = inputs(self.is_rna, sequence, structure)?;
        terms_kcal(self.model.energy_terms(&sequence, &pairings))
    }
}

/// The energy model of ViennaRNA 2.6 (dangles on both sides), without
/// pseudoknots. Default parameters: Turner 2004.
#[pyclass]
pub struct ViennaRNA {
    model: ff_energy::ViennaRNA,
    is_rna: bool,
}

#[pymethods]
impl ViennaRNA {

    #[new]
    #[pyo3(signature = (
        params = "rna_turner_2004",
        celsius = 37.0,
    ))]
    fn new(params: &str, celsius: f64) -> PyResult<Self> {
        let (thermo, is_rna) = thermo_params(params)?;
        Ok(Self { model: ff_energy::ViennaRNA::from_thermo_params(thermo, celsius), is_rna })
    }

    fn energy_of_structure(&self, sequence: &str, structure: &str) -> PyResult<f64> {
        let (sequence, pairings) = inputs(self.is_rna, sequence, structure)?;
        warn_circular(&sequence)?;
        kcal(self.model.energy_of_structure(&sequence, &pairings))
    }

    /// The energy of each loop, as `(0, description, energy in kcal/mol)`.
    /// They sum to the energy of the structure.
    fn energy_terms(&self, sequence: &str, structure: &str) -> PyResult<Vec<(usize, String, f64)>> {
        let (sequence, pairings) = inputs(self.is_rna, sequence, structure)?;
        warn_circular(&sequence)?;
        terms_kcal(self.model.energy_terms(&sequence, &pairings))
    }
}

/// The ViennaRNA model scores the open circle as 0, but charges the closure
/// of every loop of a structured circle, relative to the open chain.
fn warn_circular(sequence: &NucleotideVec) -> PyResult<()> {
    if sequence.has_circular_strand() {
        Python::attach(|py| PyErr::warn(py, &py.get_type::<pyo3::exceptions::PyUserWarning>(),
            c"ViennaRNA scores the open circle as 0, which is not the reference state of the structured circles (relative to the open chain); use Fuzzyfold.", 1))?;
    }
    Ok(())
}

/// A parameter set by name, and whether it is for RNA.
pub(crate) fn thermo_params(params: &str) -> PyResult<(&'static RNAThermoParams, bool)> {
    match params {
        "rna_turner_2004" => Ok((&RNA_TURNER_2004, true)),
        "rna_extended" => Ok((&RNA_EXTENDED, true)),
        "dna" => Ok((&DNA_MATHEWS_2004, false)),
        _ => Err(PyValueError::new_err(format!(
            "Unknown parameter set '{}'. Valid options are: 'rna_turner_2004', 'rna_extended', 'dna'.",
            params
        ))),
    }
}

/// The parsed sequence and structure, checked against each other.
fn inputs(is_rna: bool, sequence: &str, structure: &str) -> PyResult<(NucleotideVec, MultiPairTable)> {
    let sequence = match is_rna {
        true => NucleotideVec::try_from_rna(sequence),
        false => NucleotideVec::try_from_dna(sequence),
    }.map_err(|e| PyValueError::new_err(e.to_string()))?;
    let db = DotBracketVec::try_from(structure)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    check_sequence(&sequence, &db).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let pairings = MultiPairTable::try_from(&db)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok((sequence, pairings))
}

fn evaluation_error(e: EnergyError) -> PyErr {
    PyValueError::new_err(format!("Energy evaluation errored: {}", e))
}

fn kcal(energy: Result<i32, EnergyError>) -> PyResult<f64> {
    energy.map(|e| e as f64 / 100.0).map_err(evaluation_error)
}

fn terms_kcal(terms: Result<Vec<(usize, String, i32)>, EnergyError>) -> PyResult<Vec<(usize, String, f64)>> {
    let terms = terms.map_err(evaluation_error)?;
    Ok(terms.into_iter().map(|(l, t, e)| (l, t, e as f64 / 100.0)).collect())
}

fn structure_of(structure: &str) -> PyResult<MultiPairTable> {
    MultiPairTable::try_from(structure).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// A stem: its pairs, and whether it is strained.
type Stem = (Vec<(usize, usize)>, bool);

/// The stems of the pseudoknots of a structure, as the energy evaluation 
/// sees them: `(pairs, strained)` per stem, pairs outer -> inner, sorted by
/// the 5' index of the outer pair. A stem continues through stacks and
/// interior loops; it is strained if its pairs lie on a cycle of the loop
/// graph, relaxed if it crosses another pair but does not. Indices count the
/// strand breaks ('+') as positions, like the structure string. Pairs that
/// cross no pair and lie on no cycle are not listed.
#[pyfunction]
pub fn stems(structure: &str) -> PyResult<Vec<Stem>> {
    let mpt = structure_of(structure)?;
    let mut inner_of = std::collections::HashMap::new();
    let mut strained = std::collections::HashSet::new();
    // A closed loop of two pairs links them: (p, q) inside (i, j).
    let mut link = |l: &ff_energy::NestedLoop| if let (Gap::None, [(p, q), (j, i)]) = (l.gap(), l.pairs()) {
        inner_of.insert((*i, *j), (*p, *q));
    };
    for r in mpt.rings() {
        match &r {
            Ring::Loop(l) => link(l),
            Ring::Pseudoknot(pk) => {
                for l in pk.stacks.iter().chain(&pk.interior_loops) {
                    link(l);
                    strained.extend(l.pairs().iter().map(|&(a, b)| (a.min(b), a.max(b))));
                }
                for l in &pk.loops {
                    for (&(a, b), &end) in l.nested.pairs().iter().zip(&l.stem_ends) {
                        if end { strained.insert((a.min(b), a.max(b))); }
                    }
                }
            }
        }
    }
    let all: Vec<(usize, usize)> = (0..mpt.len())
        .filter_map(|i| match mpt.get(i) {
            ff_structure::MultiStruct::Paired(j) if i < j => Some((i, j)),
            _ => None,
        })
        .collect();
    let crosses = |&(i, j): &(usize, usize)| all.iter().any(|&(k, l)| i < k && k < j && j < l || k < i && i < l && l < j);
    let inner: std::collections::HashSet<_> = inner_of.values().copied().collect();
    let mut result = Vec::new();
    for &p in all.iter().filter(|p| !inner.contains(p)) {
        let mut stem = vec![p];
        while let Some(&q) = inner_of.get(stem.last().unwrap()) {
            stem.push(q);
        }
        let on_cycle = stem.iter().all(|p| strained.contains(p));
        if on_cycle || stem.iter().any(crosses) {
            result.push((stem, on_cycle));
        }
    }
    Ok(result)
}


