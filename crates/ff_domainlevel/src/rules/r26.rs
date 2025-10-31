
use crate::complexregistry::{ComplexRef, ComplexRegistry, ComplexRegistryError};
use crate::rules::RewriteRule;
use crate::complex::get_kernel;
use crate::domain::is_complement;
use structure::DotBracket;

pub struct R26;

impl RewriteRule for R26 {
    fn id(&self) -> &'static str {
        "R2.6"
    }

    fn pattern(&self) -> &'static str {
        "X( ? X ? Y) → X ? X( ? Y)"
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
            new_struc[i] = DotBracket::Unpaired;
            new_struc[j] = DotBracket::Open;
            new_struc[k] = DotBracket::Close;
            let kernel = get_kernel(&seq, &new_struc);
            ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (0..seq.len() - 2)
            .filter(move |&i| struc[i] == DotBracket::Open)
            .flat_map(move |i| {
                let k = table[i].expect("tmp"); // 0-based match
                (i + 1..k)
                    .filter(move |&j| {
                        is_complement(&seq[j], &seq[k]) &&
                            table.is_well_formed(i+1, j) &&
                            table.is_well_formed(j+1, k)
                    })
                    .filter_map(move |j| {
                        match apply_match(complex, i, j, k) {
                            Ok(product) => {
                                let rewrite = format!(
                                    "{}{}( {}{} {}) -> {}{} {}{}( {})", 
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
    fn test_enumerate_r26() {
        let i1 = ComplexRegistry::get_or_create("r26 a( b a* a )", Some("R26_I1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r26 a b a* a( )", Some("R26_O1")).expect("must be valid.");

        let my_clxs: Vec<_> = vec![o1].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a( 5a 6) -> 2a 5a( 6)"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R26.apply(&i1)
                           .map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite))
                           .unzip();

        for (product, rewrite) in R26.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

