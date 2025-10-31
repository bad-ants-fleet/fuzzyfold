use ff_structure::DotBracket;

use crate::ComplexRef;
use crate::ComplexRegistry;
use crate::ComplexRegistryError;
use crate::is_complement;
use crate::get_kernel;
use crate::rules::RewriteRule;

pub struct R23;

impl RewriteRule for R23 {
    fn id(&self) -> &'static str {
        "R2.3"
    }

    fn pattern(&self) -> &'static str {
        "X( ? Y) ? Y → X( ? Y ? Y)"
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
            let kernel = get_kernel(seq, &new_struc);
            ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (0..seq.len() - 2)
            .filter(move |&i| struc[i] == DotBracket::Open)
            .flat_map(move |i| {
                let j = table[i].expect("") as usize; // 0-based match
                (j + 1..seq.len())
                    .filter(move |&k| {
                        is_complement(&seq[i], &seq[k]) &&
                            table.is_well_formed(j+1, k)
                    })
                    .filter_map(move |k| {
                        match apply_match(complex, i, j, k) {
                            Ok(product) => {
                                let rewrite = format!(
                                    "{}{}( {}) {}{} -> {}{}( {}{} {})", 
                                    i + 1, seq[i], j + 1, k + 1, seq[k], 
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
    fn test_enumerate_r23() {
        let i1 = ComplexRegistry::get_or_create("r23 a( b ) a a*", Some("R23_I1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r23 a( b a* a )", Some("R23_O1")).expect("must be valid.");

        let my_clxs: Vec<_> = [o1].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a( 4) 6a* -> 2a( 4a* 6)"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R23.apply(&i1)
                           .map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite))
                           .unzip();

        for (product, rewrite) in R23.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

