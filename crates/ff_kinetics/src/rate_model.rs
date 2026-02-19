use crate::Move;

pub const K0: f64 = 273.15;
pub const KB: f64 = 0.001987204285; // kcal/(mol*K)

pub trait RateModel {
    /// Given dE (in kcal/mol), return the rate constant.
    fn rate(&self, m: &Move, delta_e: i32) -> f64;
}

#[derive(Debug, Clone, Copy)]
pub struct Metropolis {
    kt: f64, // k_B * T in kcal/mol
    k0: f64,
    k3ws: f64,
    k4ws: f64,
}

impl Metropolis {
    pub fn new(celsius: f64, k0: f64, k3ws: Option<f64>, k4ws: Option<f64>) -> Self {
        if k0 < 0. {
            panic!("k0 must not be negative!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            kt: KB * t_kelvin,
            k0,
            k3ws: k3ws.unwrap_or(0.0),
            k4ws: k4ws.unwrap_or(0.0),
        }
    }
}

impl RateModel for Metropolis {
    fn rate(&self, mv: &Move, delta_e: i32) -> f64 {
        match &mv {
            Move::Add { .. } | Move::Del { .. } => {
                if delta_e <= 0 {
                    self.k0
                } else {
                    self.k0 * ((-delta_e as f64 / 100.) / self.kt).exp()
                }
            },
            Move::ShiftIK { .. } | Move::ShiftJK { .. } => {
                if delta_e <= 0 {
                    self.k3ws
                } else {
                    self.k3ws * ((-delta_e as f64 / 100.) / self.kt).exp()
                }
            },
            Move::ShiftIKLJ { .. } | Move::ShiftILJK { .. } => {
                if delta_e <= 0 {
                    self.k4ws
                } else {
                    self.k4ws * ((-delta_e as f64 / 100.) / self.kt).exp()
                }
            },
        }
   }
}    

#[derive(Debug, Clone, Copy)]
pub struct Kawasaki {
    beta: f64, // 2 * k_B * T in kcal/mol
    k0: f64,
}

impl Kawasaki {
    pub fn new(celsius: f64, k0: f64) -> Self {
        if k0 < 0. {
            panic!("k0 must not be negative!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            beta: KB * t_kelvin * 2.0,
            k0,
        }
    }
}

impl RateModel for Kawasaki {
    fn rate(&self, mv: &Move, delta_e: i32) -> f64 {
        match &mv {
            Move::Add { .. } | Move::Del { .. } => {
                self.k0 * ((-delta_e as f64 / 100.) / self.beta).exp()
            },
            _ => unreachable!("Kawasaki does not support shift activation energies.")
        }

    }
}

