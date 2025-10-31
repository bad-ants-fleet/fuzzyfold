
use crate::complexregistry::{ComplexRef, ComplexRegistry, ComplexRegistryError};
use crate::rules::RewriteRule;
use crate::domain::is_complement;
use crate::complex::get_kernel;
use structure::DotBracket;

pub struct R25;

impl RewriteRule for R25 {
    fn id(&self) -> &'static str {
        "R2.5"
    }

    fn pattern(&self) -> &'static str {
        "X ? X( ? Y) → X( ? X ? Y)"
    }

    fn category(&self) -> &'static str {
        "three-way branch migration"
    }

    fn cost(&self) -> f64 {
        1.2
    }

    fn apply<'a>(&'a self, complex: &'a ComplexRef) -> 
        Box<dyn Iterator<Item = (ComplexRef, String)> + 'a> {
        let seq = complex.sequence();
        let struc = complex.structure();
        let table = complex.pair_table();

        fn apply_match(complex: &ComplexRef, i:usize, j:usize, k:usize) 
            -> Result<ComplexRef, ComplexRegistryError> {
            let seq = complex.sequence();
            let mut new_struc = complex.structure().to_vec();
            new_struc[i] = DotBracket::Open;
            new_struc[j] = DotBracket::Unpaired;
            new_struc[k] = DotBracket::Close;
            let kernel = get_kernel(&seq, &new_struc);
            ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (1..seq.len()).rev()
            .filter(move |&k| struc[k] == DotBracket::Close)
            .flat_map(move |k| {
                let j = table[k].expect(""); // 0-based match
                (0..j)
                    .filter(move |&i| {
                        is_complement(&seq[i], &seq[k]) &&
                            table.is_well_formed(i+1, j)
                    })
                    .filter_map(move |i| {
                        match apply_match(complex, i, j, k) {
                            Ok(product) => {
                                let rewrite = format!(
                                    "{}{} {}{}( {}) -> {}{}( {}{} {})", 
                                    i + 1, seq[i], j + 1, seq[j], k + 1, 
                                    i + 1, seq[i], j + 1, seq[j], k + 1);
                                Some((product, rewrite))
                            }
                            Err(_) => None,
                        }
                    })
            })
        )
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_enumerate_r25() {
        let i1 = ComplexRegistry::get_or_create("r25 a b a* a( )", Some("R25_I1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r25 a( b a* a )", Some("R25_O1")).expect("must be valid.");

        let my_clxs: Vec<_> = vec![o1].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a 5a( 6) -> 2a( 5a 6)"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R25.apply(&i1)
                           .map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite))
                           .unzip();

        for (product, rewrite) in R25.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

