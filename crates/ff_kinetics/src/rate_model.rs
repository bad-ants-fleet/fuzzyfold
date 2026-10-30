use crate::Move;
use crate::MoveTypes;

/// 0 degrees Celsius in Kelvin.
pub const ZERO_CELSIUS: f64 = 273.15;
/// The gas constant in kcal/(mol K).
pub const KB: f64 = 0.001987204285;

pub trait RateModel: Clone + Send + Sync {
    /// Given dE (in kcal/mol), return the rate constant.
    fn rate(&self, m: &Move, delta_e: i32) -> f64;
}

/// The Arrhenius rate model.
///
/// One rate constant per move type: k0 for adding and deleting pairs, k3ws
/// for three-way and k4ws for four-way shift moves, and k_atoi for A-to-I
/// edits (independent of the energy: editing is irreversible). A rate
/// constant of 0 turns that move type off (see `move_types`: the move set
/// then need not propose them).
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
    /// The rate of an A-to-I edit at a site.
    k_atoi: f64,
}

impl Arrhenius {
    pub fn new(celsius: f64, k0: f64, k3ws: Option<f64>, k4ws: Option<f64>, k_atoi: Option<f64>) -> Self {
        if k0 < 0. {
            panic!("k0 must not be negative!");
        }
        let t_kelvin = celsius + ZERO_CELSIUS;
        Self { 
            kt: KB * t_kelvin,
            k0,
            k3ws: k3ws.unwrap_or(0.0),
            k4ws: k4ws.unwrap_or(0.0),
            k_atoi: k_atoi.unwrap_or(0.0),
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

    pub fn k_atoi(&self) -> Option<f64> {
        (self.k_atoi > 0.0).then_some(self.k_atoi)
    }

    /// The move types with a rate: the move set to simulate with.
    pub fn move_types(&self) -> MoveTypes {
        MoveTypes { three_way: self.k3ws > 0.0, four_way: self.k4ws > 0.0, a_to_i: self.k_atoi > 0.0 }
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
            Move::Edit { .. } => self.k_atoi,
        }
   }
}    

