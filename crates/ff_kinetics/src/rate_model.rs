use crate::Move;

pub const K0: f64 = 273.15;
pub const KB: f64 = 0.001987204285; // kcal/(mol*K)

pub trait RateModel: Send + Sync {
    /// Given dE (in kcal/mol), return the rate constant.
    fn rate(&self, m: &Move, delta_e: i32) -> f64;
}

/// The Arrhenius rate model.
///
/// We specify kT = k_Boltzmann * temperature [kcal/mol].
/// There are multiple types of k0 available to parameterize different move
/// sets. At the moment, we support the typcal add/delete moves: k0, three-way
/// shit moves: k3ws and four-way shit moves: k4ws. Setting either of those to 0
/// effictively turns the move set off. (However, for performancy reasons you
/// may want to switch these moves sets off during neighborhood generation as
/// well.)
#[derive(Debug, Clone, Copy)]
pub struct Arrhenius {
    /// kT = k_Boltzmann * temperature [kcal/mol].
    kt: f64,
    /// The maximum rate for base-pair formation/breaking (k_0 = A_0 * exp(-G_0/kT)) 
    k0: f64,
    /// The maximum rate for three-way shift moves (k_{3ws} = A_{3ws} * exp(-G_{3ws}/kT)) 
    k3ws: f64,
    /// The maximum rate for four-way shift moves (k_{4ws} = A_{4ws} * exp(-G_{4ws}/kT)) 
    k4ws: f64,
}

impl Arrhenius {
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

    pub fn k0(&self) -> Option<f64> {
        if self.k0 > 0.0 {
            Some(self.k0)
        } else {
            None
        }
    }

    pub fn k3ws(&self) -> Option<f64> {
        if self.k3ws > 0.0 {
            Some(self.k3ws)
        } else {
            None
        }
    }

    pub fn k4ws(&self) -> Option<f64> {
        if self.k4ws > 0.0 {
            Some(self.k4ws)
        } else {
            None
        }
    }

}

impl RateModel for Arrhenius {
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

