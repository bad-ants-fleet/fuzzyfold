
use crate::NearestNeighborLoop;
use crate::LoopDecomposition;
use crate::Base;

pub const K0: f64 = 273.15;

pub trait EnergyModel {
    fn can_pair(&self, b1: Base, b2: Base) -> bool;

    fn min_hairpin_size(&self) -> usize;

    fn temperature(&self) -> f64;

    fn energy_of_structure<T: LoopDecomposition>(&self, 
        sequence: &[Base], 
        structure: &T
    ) -> i32;

    //TODO: consider returning Result in case the evaluation fails!
    fn energy_of_loop(&self, 
        sequence: &[Base], 
        nn_loop: &NearestNeighborLoop
    ) -> i32;
}

#[cfg(test)]
mod tests {
    use super::*;
    use Base::*;

    struct MockEnergyModel;

    impl EnergyModel for MockEnergyModel {
        fn can_pair(&self, b1: Base, b2: Base) -> bool {
            matches!((b1, b2), (A, U) | (U, A) | (C, G) | (G, C))
        }
        
        fn min_hairpin_size(&self) -> usize {
            3
        }

        fn temperature(&self) -> f64 {
            37.0
        }

        fn energy_of_structure<T: LoopDecomposition>(
            &self,
            _sequence: &[Base],
            _structure: &T,
        ) -> i32 {
            -10
        }

        fn energy_of_loop(
            &self,
            _sequence: &[Base],
            _nn_loop: &NearestNeighborLoop,
        ) -> i32 {
            5 
        }
    }

    #[test]
    fn test_can_pair() {
        let model = MockEnergyModel;
        assert!(model.can_pair(A, U));
        assert!(model.can_pair(C, G));
        assert!(!model.can_pair(A, G));
        assert!(!model.can_pair(C, C));
    }

    #[test]
    fn test_min_hairpin_size() {
        let model = MockEnergyModel;
        assert_eq!(model.min_hairpin_size(), 3);
    }

    #[test]
    fn test_energy_of_loop() {
        let model = MockEnergyModel;

        let sequence = vec![A, U, C, G];
        let nn_loop = NearestNeighborLoop::Hairpin { closing: (0, 3) };

        let energy = model.energy_of_loop(&sequence, &nn_loop);
        assert_eq!(energy, 5);
    }
}


