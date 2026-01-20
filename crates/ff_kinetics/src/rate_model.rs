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
}

impl Metropolis {
    pub fn new(celsius: f64, k0: f64) -> Self {
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            kt: KB * t_kelvin,
            k0,
        }
    }
}

impl RateModel for Metropolis {
    fn rate(&self, _: &Move, delta_e: i32) -> f64 {
        if delta_e <= 0 {
            self.k0
        } else {
            self.k0 * ((-delta_e as f64 / 100.) / self.kt).exp()
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
        if k0 <= 0. {
            panic!("k0 must be positive!");
        }
        let t_kelvin = celsius + K0;
        Self { 
            beta: KB * t_kelvin * 2.0,
            k0,
        }
    }
}

impl RateModel for Kawasaki {
    fn rate(&self, _: &Move, delta_e: i32) -> f64 {
        self.k0 * ((-delta_e as f64 / 100.) / self.beta).exp()
    }
}

//pub struct MixedModel {
//    kt: f64,      // k_B * T in kcal/mol 
//    k_m: f64,     // Metropolis rate constant
//    k_k3way: f64, // Kawasaki three-way shift
//    k_k4way: f64, // Kawasaki four-way shift
//}
//
//impl MixedModel {
//    pub fn new(celsius: f64, k_m: f64, k_k3way: f64, k_k4way: f64) -> Self {
//        let t_kelvin = celsius + K0;
//        Self { 
//            kt: KB * t_kelvin,
//            k_m,
//            k_k3way,
//            k_k4way,
//        }
//    }
//
//    fn mrate(&self, delta_e: i32) -> f64 {
//        if delta_e <= 0 {
//            self.k_m
//        } else {
//            self.k_m * ((-delta_e as f64 / 100.) / self.kt).exp()
//        }
//    }
//}
//
//impl RateModel for MixedModel{
//    fn rate(&self, mv: &Move, delta_e: i32) -> f64 {
//        match mv {
//            _ => self.mrate(delta_e),
//        }
//    }
//}

