use crate::parameters::LoopEntry;

pub trait EnergyParameter: Sized {
    type Output;
    fn rescale(&self, enthalpy: &Self, scale: f64) -> Self::Output;
}

impl EnergyParameter for i32 {
    type Output = i32;
    #[inline]
    fn rescale(&self, h: &Self, scale: f64) -> Self {
        let g37 = *self as f64;
        let h = *h as f64;
        (h - (h - g37) * scale).round() as i32
    }
}

impl<T: EnergyParameter<Output = T> + Copy, const N: usize> EnergyParameter for [T; N] {
    type Output = [T; N];
    #[inline]
    fn rescale(&self, enthalpies: &Self, scale: f64) -> Self {
        let mut new = *self;
        for (g, h) in new.iter_mut().zip(enthalpies.iter()) {
            *g = g.rescale(h, scale);
        }
        new
    }
}

impl EnergyParameter for LoopEntry {
    type Output = LoopEntry;
    #[inline]
    fn rescale(&self, enthalpy: &Self, scale: f64) -> Self {
        Self {
            seq: self.seq,
            val: self.val.rescale(&enthalpy.val, scale),
        }
    }
}

impl<T: EnergyParameter<Output = T> + Clone> EnergyParameter for &[T] {
    type Output = Vec<T>;
    #[inline]
    fn rescale(&self, enthalpies: &Self, scale: f64) -> Vec<T> {
        self.iter()
            .zip(enthalpies.iter())
            .map(|(g, h)| g.rescale(h, scale))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const T_REF: f64 = 310.15;
    use crate::K0;

    #[test]
    fn rescale_i32() {
        let g37 = -210;
        let h = -1050;
        let scale = T_REF / T_REF; // = 1.0
        let gt = g37.rescale(&h, scale);
        assert_eq!(gt, g37);
    }

    #[test]
    fn rescale_i32_known_value() {
        let g37 = -210;
        let h = -1050;

        // Example: 25°C
        let kelvin = 25.0 + K0;
        let scale = kelvin / T_REF;
        let gt = g37.rescale(&h, scale);

        // Manual computation
        let expected = {
            let g = g37 as f64;
            let h = h as f64;
            (h - (h - g) * scale).round() as i32
        };
        assert_eq!(gt, expected);
    }

    #[test]
    fn rescale_array() {
        let g37 = [[-210, -220], [-330, -340]];
        let h = [[-1050, -1140], [-1340, -1490]];

        let kelvin = 25.0 + K0;
        let scale = kelvin / T_REF;

        let gt = g37.rescale(&h, scale);

        for i in 0..2 {
            for j in 0..2 {
                let expected =
                    (h[i][j] as f64
                        - (h[i][j] as f64 - g37[i][j] as f64) * scale)
                        .round() as i32;
                assert_eq!(gt[i][j], expected);
            }
        }
    }

}
