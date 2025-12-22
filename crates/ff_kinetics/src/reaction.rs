use ff_structure::NAIDX;
use crate::RateModel;

#[derive(Debug, Clone, PartialEq)]
pub enum Reaction {
    Add {
        i: NAIDX,
        j: NAIDX,
        delta_e: i32,
        log_rate: f64,
    },
    Del {
        i: NAIDX,
        j: NAIDX,
        delta_e: i32,
        log_rate: f64,
    },
}

impl Reaction {
    pub fn new_add<K: RateModel>(model: &K, 
        i: NAIDX, j: NAIDX, delta_e: i32
) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Add { i, j, delta_e, log_rate: rate }
    }

    pub fn new_del<K: RateModel>(model: &K, 
        i: NAIDX, j: NAIDX, delta_e: i32) -> Self {
        let rate = model.log_rate(delta_e);
        Reaction::Del { i, j, delta_e, log_rate: rate }
    }

    pub fn ij(&self) -> (NAIDX, NAIDX) {
        match self {
            Reaction::Add { i, j, .. } => (*i, *j),
            Reaction::Del { i, j, .. } => (*i, *j),
        }
    }

    pub fn log_rate(&self) -> f64 {
        match self {
            Reaction::Add { log_rate, .. } => *log_rate,
            Reaction::Del { log_rate, .. } => *log_rate,
        }
    }

    pub fn delta_e(&self) -> i32 {
        match self {
            Reaction::Add { delta_e, .. } => *delta_e,
            Reaction::Del { delta_e, .. } => *delta_e,
        }
    }

}

