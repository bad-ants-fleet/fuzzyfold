use ff_structure::DotBracket;

use crate::ComplexRef;
use crate::ComplexRegistry;
use crate::ComplexRegistryError;
use crate::is_complement;
use crate::get_kernel;
use crate::rules::RewriteRule;

pub struct R12;

impl RewriteRule for R12 {
    fn id(&self) -> &'static str {
        "R1.2"
    }

    fn pattern(&self) -> &'static str {
        "X( ? ) => X ? Y"
    }

    fn category(&self) -> &'static str {
        "unbind"
    }

    fn cost(&self) -> f64 {
        1.5  // Arbitrary example, tune as needed
    }

    fn apply<'a>(&'a self, complex: &'a ComplexRef) -> 
        Box<dyn Iterator<Item = (ComplexRef, String)> + 'a> {
        let seq = complex.sequence();
        let struc = complex.structure();
        let table = complex.pair_table();

        fn apply_match(complex: &ComplexRef, i: usize, j: usize) -> 
            Result<ComplexRef, ComplexRegistryError> {
                let seq = complex.sequence();
                let mut new_struc = complex.structure().to_vec();
                new_struc[i] = DotBracket::Unpaired;
                new_struc[j] = DotBracket::Unpaired;
                let kernel = get_kernel(seq, &new_struc);
                ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (0..seq.len()-1)
            .filter(move |&i| struc[i] == DotBracket::Open)
            .filter_map(move |i| {
                let j = table[i].expect("") as usize;
                assert!(j > i, "in r12: j < i ({} < {})", j, i);
                assert!(seq.len() > j, "in r12: j > len ({} > {})", j, seq.len());
                assert_eq!(struc[j], DotBracket::Close, "({:?} > {:?})", struc[j], DotBracket::Close);
                assert!(is_complement(&seq[i], &seq[j]));
                match apply_match(complex, i, j) {
                    Ok(product) => {
                        let rewrite = format!(
                            "{}{}( {}) -> {}{} {}{}", 
                            i + 1, seq[i], j + 1,
                            i + 1, seq[i], j + 1, seq[j]);
                        Some((product, rewrite))
                    }
                    Err(_) => None,
                }
            })
        )
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enumerate_r12() {
        let i1 = ComplexRegistry::get_or_create("r12 a( b ) a( )", Some("R12_I1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r12 a b a* a( )", Some("R12_O2")).expect("must be valid.");
        let o2 = ComplexRegistry::get_or_create("r12 a( b ) a a*", Some("R12_O1")).expect("must be valid.");

        let my_clxs: Vec<_> = [o1, o2].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a( 4) -> 2a 4a*", "5a( 6) -> 5a 6a*"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R12.apply(&i1)
                           .map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite))
                           .unzip();

        for (product, rewrite) in R12.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

