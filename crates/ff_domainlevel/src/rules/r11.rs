use ff_structure::DotBracket;

use crate::ComplexRef;
use crate::ComplexRegistry;
use crate::ComplexRegistryError;
use crate::is_complement;
use crate::get_kernel;
use crate::rules::RewriteRule;

pub struct R11;

impl RewriteRule for R11 {
    fn id(&self) -> &'static str {
        "R1.1"
    }

    fn pattern(&self) -> &'static str {
        "X ? Y → X( ? )"
    }

    fn category(&self) -> &'static str {
        "bind"
    }

    fn cost(&self) -> f64 {
        1.5 // arbitrary
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
                new_struc[i] = DotBracket::Open;
                new_struc[j] = DotBracket::Close;
                let kernel = get_kernel(seq, &new_struc);
                ComplexRegistry::get_or_create(&kernel, None)
        }

        Box::new(
            (0..seq.len()-1)
            .filter(move |&i| struc[i] == DotBracket::Unpaired)
            .flat_map(move |i| {
                (i + 1..seq.len())
                    .filter(move |&j| {
                        struc[j] == DotBracket::Unpaired &&
                            is_complement(&seq[i], &seq[j]) &&
                            table.is_well_formed(i+1, j)
                    })
                    .filter_map(move |j| {
                        match apply_match(complex, i, j) {
                            Ok(product) => {
                                let rewrite = format!(
                                    "{}{} {}{} -> {}{}( {})", 
                                    i + 1, seq[i], j + 1, seq[j],
                                    i + 1, seq[i], j + 1);
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
    fn test_enumerate_r11() {
        let i1 = ComplexRegistry::get_or_create("r11 a b a* a a*", Some("TI1")).expect("must be valid.");
        let o1 = ComplexRegistry::get_or_create("r11 a( b ) a a*", Some("TO1")).expect("must be valid.");
        let o2 = ComplexRegistry::get_or_create("r11 a( b a* a )", Some("TO2")).expect("must be valid.");
        let o3 = ComplexRegistry::get_or_create("r11 a b a*( ) a*", Some("TO3")).expect("must be valid.");
        let o4 = ComplexRegistry::get_or_create("r11 a b a* a( )", Some("TO4")).expect("must be valid.");

        let my_clxs: Vec<_> = [o1, o2, o3, o4].iter().map(|c| c.kernel().to_string()).collect();
        let my_rwrs = vec!["2a 4a* -> 2a( 4)", "2a 6a* -> 2a( 6)", "4a* 5a -> 4a*( 5)", "5a 6a* -> 5a( 6)"];

        let (clxs, rwrs): (Vec<_>, Vec<_>) = R11.apply(&i1).map(|(cplx, rewrite)| (cplx.kernel().to_string(), rewrite)).unzip();

        for (product, rewrite) in R11.apply(&i1) {
            println!("{} -> {} ({})", i1, product, rewrite);
        }
        assert_eq!(my_clxs, clxs);
        assert_eq!(my_rwrs, rwrs);
    }
} 

