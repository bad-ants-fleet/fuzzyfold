use std::collections::{HashSet, VecDeque};

use crate::complexregistry::ComplexRef;
use crate::reactions::Reaction;
use crate::rules::RewriteRule;
use crate::rules::*;

pub fn enumerate_complexes(
    start: impl IntoIterator<Item = ComplexRef>,
    rules: &[&dyn RewriteRule],
) -> (Vec<ComplexRef>, Vec<Reaction>) {
    let mut visited = HashSet::new();
    let mut queue = VecDeque::new();
    let mut complexes = Vec::new();
    let mut reactions = Vec::new();

    for complex in start {
        let k = complex.kernel().to_string();
        if visited.insert(k) {
            queue.push_back(complex.clone());
            complexes.push(complex);
        }
    }

    while let Some(current) = queue.pop_front() {
        for rule in rules {
            for reaction in rule.get_reactions(&current) {
                let product = reaction.product();
                let key = product.kernel().to_string();
                if visited.insert(key) {
                    queue.push_back(product.clone());
                    complexes.push(product.clone());
                }
                reactions.push(reaction);
            }
        }
    }

    (complexes, reactions)
}

pub fn enumerate_all_rules(
    start: impl IntoIterator<Item = ComplexRef>,
) -> (Vec<ComplexRef>, Vec<Reaction>) {
    enumerate_complexes(start, all_rules())
} 

pub fn enumerate_bind_threeway_fourway(
    start: impl IntoIterator<Item = ComplexRef>,
) -> (Vec<ComplexRef>, Vec<Reaction>) {
    enumerate_complexes(start, bind_threeway_fourway())
} 

pub fn enumerate_bind_threeway(
    start: impl IntoIterator<Item = ComplexRef>,
) -> (Vec<ComplexRef>, Vec<Reaction>) {
    enumerate_complexes(start, bind_threeway())
} 


#[cfg(test)]
mod tests {
    use super::*;
    use crate::ComplexRegistry;
    use crate::rules::r11::R11;

    #[test]
    fn test_enumerate_single_rule() {
        let i1 = ComplexRegistry::get_or_create("en a b a* a a*", Some("EN_I1")).expect("must be valid.");
        
        let start = vec![i1.clone()];
        let rules: &[&dyn RewriteRule] = &[&R11]; // etc.

        let (cxs, rxs) = enumerate_complexes(start.clone(), rules);

        println!("Enumerated {} complexes starting from '{:?}'", cxs.len(), start.iter().map(|c| c.kernel().to_string()).collect::<Vec<_>>());

        for cx in &cxs {
            println!("{}", cx);
        }

        assert_eq!(cxs.len(), 7);
        assert_eq!(rxs.len(), 8);
    }

    #[test]
    fn test_enumerate_all_rules() {
        let i1 = ComplexRegistry::get_or_create("enar a b a* a a*", Some("ENAR_I1")).expect("must be valid.");
        
        let start = vec![i1.clone()];
        let (cxs, rxs) = enumerate_all_rules(start.clone());

        println!("Enumerated {} complexes starting from '{:?}'",
            cxs.len(), start.iter().map(|c| c.kernel().to_string()).collect::<Vec<_>>());

        for rx in &rxs {
            println!("{}", rx);
        }

        assert_eq!(cxs.len(), 7);
        assert_eq!(rxs.len(), 26); // not exhaustively verified
    }
} 

