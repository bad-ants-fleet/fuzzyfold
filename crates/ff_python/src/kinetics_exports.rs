use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::sync::Arc;
use rand::SeedableRng;
use rand::rngs::SmallRng;

use ff_structure::DotBracketVec;
use ff_structure::PairTable;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA;
use ff_kinetics::SSA;
use ff_kinetics::shift_policy;
use ff_kinetics::Arrhenius;
use ff_kinetics::Walker;
use ff_kinetics::LoopNeighbors;

//TODO: support shifts, rename to arrhenius

#[pyclass]
pub struct Simulator {
    energy_model: Arc<ViennaRNA>,
    rate_model: Arrhenius,
}


#[pymethods]
impl Simulator {
    #[new]
    #[pyo3(signature = (
        celsius=37.0,
        k0=1e5,
        k3ws=0.0,
        k4ws=0.0,
    ))]
    fn new(
        celsius: f64,
        k0: f64,
        k3ws: f64,
        k4ws: f64,
    ) -> PyResult<Self> {

        if k0 < 0.0 || k3ws < 0.0 || k4ws < 0.0 {
            return Err(PyValueError::new_err(
                "Rate constants must be non-negative",
            ));
        }

        let mut energy_model = ViennaRNA::default();
        energy_model.reset_with_temperature(celsius);

        let rate_model = Arrhenius::new(
            celsius,
            k0,
            Some(k3ws),
            Some(k4ws),
        );

        Ok(Self {
            energy_model: Arc::new(energy_model),
            rate_model,
        })
    }

    #[pyo3(signature = (
            sequence,
            start=None,
            t_ext=None,
            t_end=1.0,
    ))]
    fn simulate(
        &self,
        sequence: &str,
        start: Option<&str>,
        t_ext: Option<f64>,
        t_end: f64,
    ) -> PyResult<SimulationIterator> {

        let seq = NucleotideVec::try_from(sequence)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;

        let start_db = match start {
            Some(s) => DotBracketVec::try_from(s)
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
            None => DotBracketVec::try_from(".")
                .map_err(|e| PyValueError::new_err(e.to_string()))?,
        };

        if start_db.len() < seq.len() && t_ext.is_none() {
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

        match (self.rate_model.k3ws().is_some(), self.rate_model.k4ws().is_some()) {
            (false, false) => build_iterator(
                seq,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::NoShift,
                SSAKind::NoShift,
            ),

            (true, false) => build_iterator(
                seq,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeWayOnly,
                SSAKind::ThreeWayOnly,
            ),

            (false, true) => build_iterator(
                seq,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::FourWayOnly,
                SSAKind::FourWayOnly,
            ),

            (true, true) => build_iterator(
                seq,
                &start_pt,
                Arc::clone(&self.energy_model),
                self.rate_model,
                times,
                shift_policy::ThreeAndFour,
                SSAKind::ThreeAndFour,
            ),
        }
   }
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
                            // advance the simulator to update the structure.
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


