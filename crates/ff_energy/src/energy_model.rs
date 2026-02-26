use std::fmt;

use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;

pub const K0: f64 = 273.15;

#[derive(Debug)]
pub enum EnergyError {
    HairpinTooSmall { size: usize, min: usize },
    InvalidClosingPair,
}

impl fmt::Display for EnergyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            EnergyError::HairpinTooSmall { size, min } => {
                write!(f, "Hairpin too small: size {}, minimum allowed is {}", size, min)
            }
            EnergyError::InvalidClosingPair => {
                write!(f, "Invalid closing base pair")
            }
        }
    }
}

impl std::error::Error for EnergyError {}

pub trait EnergyModel {
    fn can_pair(&self, b1: Base, b2: Base) -> bool;

    fn min_hairpin_size(&self) -> usize;

    fn temperature(&self) -> f64;

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> Result<i32, EnergyError>;

    fn energy_of_loop(&self,
        sequence: &[Base],
        nn_loop: &NearestNeighborLoop,
    ) -> Result<i32, EnergyError>;
}

