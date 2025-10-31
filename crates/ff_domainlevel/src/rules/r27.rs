use ff_structure::DotBracket;

use crate::ComplexRef;
use crate::ComplexRegistry;
use crate::ComplexRegistryError;
use crate::get_kernel;
use crate::rules::RewriteRule;

pub struct R27;

impl RewriteRule for R27 {
    fn id(&self) -> &'static str {
        "R2.7"
    }

    fn pattern(&self) -> &'static str {
        "X( ? Y) ? X( ? Y) → X( ? Y( ? X) ? Y)"
    }

    fn category(&self) -> &'static str {
        "four-way branch migration"
    }

    fn cost(&self) -> f64 {
        1.2
    }

    fn apply<'a>(&'a self, complex: &'a ComplexRef) -> 
        Box<dyn Iterator<Item = (ComplexRef, String)> + 'a> {
        let seq = complex.sequence();
        let struc = complex.structure();
        let table = complex.pair_table();

        fn apply_match(complex: &ComplexRef, i:usize, j:usize, k:usize, l:usize) 
            -> Result<ComplexRef, ComplexRegistryError> {
            let seq = complex.sequence();
            let mut new_struc = complex.structure().to_vec();
            new_struc[i] = DotBracket::Open;
            new_struc[j] = DotBracket::Open;
            new_struc[k] = DotBracket::Close;
            new_struc[l] = DotBracket::Close;
            let kernel = get_kernel(seq, &new_struc);
            ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (0..seq.len()-3)
            .filter(move |&i| struc[i] == DotBracket::Open)
            .flat_map(move |i| {
                let j = table[i].expect("") as usize; // 0-based match
                (j + 1..seq.len()-1)
                    .filter(move |&k| struc[k] == DotBracket::Open)
                    .filter_map(move |k| {
                        let l = table[k].expect("") as usize; // 0-based match
                        match apply_match(complex, i, j, k, l) {
                            Ok(product) => {
                                let rewrite = format!(
                                    "{}{}( {}) {}{}( {}) -> {}{}( {}{}( {}) {})", 
                                    i + 1, seq[i], j + 1, k + 1, seq[k], l + 1,
                                    i + 1, seq[i], j + 1, seq[j], k + 1, l + 1);
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
    fn test_enumerate_r27() {
        let i1 = ComplexRegistry::get_or_create("r27 a( b ) a( )", Some("R27_I1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r27 a( b a*( ) )", Some("R27_O1")).expect("must be valid.");

        let my_clxs: Vec<_> = [o1].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a( 4) 5a( 6) -> 2a( 4a*( 5) 6)"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R27.apply(&i1)
                           .map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite))
                           .unzip();

        for (product, rewrite) in R27.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

