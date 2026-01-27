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
    ks: f64,
}

impl Metropolis {
    pub fn new(celsius: f64, k0: f64, ks: f64) -> Self {
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        if ks <= 0. {
            panic!("ks must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            kt: KB * t_kelvin,
            k0,
            ks,
        }
    }
}

impl RateModel for Metropolis {
    fn rate(&self, mv: &Move, delta_e: i32) -> f64 {
        match &mv {
            Move::ShiftIK { .. } | Move::ShiftJK { .. } => {
                if delta_e <= 0 {
                    self.ks
                } else {
                    self.ks * ((-delta_e as f64 / 100.) / self.kt).exp()
                }
            },
            Move::Add { .. } | Move::Del { .. } => {
                if delta_e <= 0 {
                    self.k0
                } else {
                    self.k0 * ((-delta_e as f64 / 100.) / self.kt).exp()
                }
            },
        }
   }
}    

#[derive(Debug, Clone, Copy)]
pub struct Kawasaki {
    beta: f64, // 2 * k_B * T in kcal/mol
    k0: f64,
    ks: f64,
}

impl Kawasaki {
    pub fn new(celsius: f64, k0: f64, ks: f64) -> Self {
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        if ks <= 0. {
            panic!("ks must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            beta: KB * t_kelvin * 2.0,
            k0,
            ks,
        }
    }
}

impl RateModel for Kawasaki {
    fn rate(&self, mv: &Move, delta_e: i32) -> f64 {
        match &mv {
            Move::ShiftIK { .. } | Move::ShiftJK { .. } => {
                self.ks * ((-delta_e as f64 / 100.) / self.beta).exp()
            },
            Move::Add { .. } | Move::Del { .. } => {
                self.k0 * ((-delta_e as f64 / 100.) / self.beta).exp()
            },
        }

    }
}

