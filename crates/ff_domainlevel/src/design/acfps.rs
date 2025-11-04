//! Abstract cotranscriptional folding pathways.


use std::fmt;
use std::error;
use std::collections::HashMap;
use ndarray::Array2;
use ahash::AHashMap;
use ahash::AHashSet;

use ff_structure::NAIDX;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;
use ff_structure::PairTable;

use crate::NussinovDP;
use crate::Pair;
use crate::PairSet;

use crate::design::{
    PartialOrder, 
    UnionFind};
use crate::DomainRegistry;
use crate::design::SegmentSequence;

#[derive(Debug, Clone)]
pub struct Acfp {
    db_path: Vec<DotBracketVec>, // That *is* the path.
    // The following fields are commonly computed / cached properties 
    // that remain valid if the path gets extended.
    pt_path: Vec<PairTable>,
    all_pairs: AHashSet<Pair>, // by convention: i < j
    partial_order: Option<PartialOrder>,
    length: usize,
}

impl Acfp {
    pub fn path(&self) -> &[DotBracketVec] {
        &self.db_path
    }

    pub fn pt_path(&self) -> &[PairTable] {
        &self.pt_path
    }

    pub fn all_pairs(&self) -> &AHashSet<Pair> {
        &self.all_pairs
    }

    pub fn pair_hierarchy(&self) -> Option<AHashMap<(NAIDX, NAIDX), usize>> {
        self.partial_order.as_ref().map(|po| po.pair_hierarchy())
    }

    pub fn extend_by_one(&mut self, new_db: DotBracketVec) {
        assert_eq!(new_db.len(), self.length + 1);
        let new_pt = PairTable::try_from(&new_db).unwrap();
        let plist = PairSet::from(&new_pt);
        for pair in plist.iter() {
            self.all_pairs.insert(pair);
        }
        
        if let Some(ref mut po) = self.partial_order {
            if !po.extend_by_pairtable(&new_pt) {
                self.partial_order = None;
            }
        }
 
        self.pt_path.push(new_pt);
        self.db_path.push(new_db);
        self.length += 1;
   }
 
    /// Every target structure is the unique best structure.
    pub fn has_valid_pair_hierarchy(&self) -> bool {
        let Some(ph) = self.pair_hierarchy() else {
            return false;
        };
       
        //NOTE: this is matix only considers explicit 
        //complementarity, which may lead to unsaturated structures.
        let n = self.length;
        let mut p = Array2::from_elem((n, n), 0);
        for (e, c) in ph {
            // We are putting 10 extra, to reward the fact 
            // that a pair is formed.
            p[(e.0 as usize, e.1 as usize)] += 100 + c;
            p[(e.1 as usize, e.0 as usize)] += 100 + c;
        }
        let ndp = NussinovDP::from(p);

        for ss in self.db_path.iter() {
            let nss = ndp.all_mfe_structs(Some(ss.len()));
            if nss.len() > 1 {
                return false
            } 
            if &nss[0] != ss {
                return false;
            }
        }
        true
    }

    pub fn all_total_orders(&self) -> Option<Vec<Vec<Pair>>> {
        self.partial_order.as_ref().map(|po| po.all_total_orders())
    }
    
    pub fn is_cycle_free(&self) -> bool {
        let mut uf = UnionFind::new(self.length + 1);
        self.all_pairs.iter().all(|&p| uf.union(p.i() as usize, p.j() as usize))
    }
 
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let mut uf = UnionFind::new(self.length);
        self.all_pairs.iter().all(|&p| uf.union(p.i() as usize, p.j() as usize));
        uf.connected_components()
    }

    pub fn has_valid_design(&self, registry: &mut DomainRegistry) -> bool {
        let Some(_) = self.pair_hierarchy() else {
            return false;
        };
        let segseq = SegmentSequence::design_from_acfp(self, registry).unwrap();
        segseq.implements_acfp(self.path(), registry)
    }

    pub fn is_valid(&self, registry: &mut DomainRegistry) -> bool {
        self.is_cycle_free() 
            && self.has_valid_pair_hierarchy()
            && self.has_valid_design(registry)
    }

}

#[derive(Debug)]
pub struct AcfpParseError {
    pub message: String,
}

impl fmt::Display for AcfpParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Error parsing ACFP: {}", self.message)
    }
}

impl error::Error for AcfpParseError {}

impl TryFrom<&str> for Acfp {
    type Error = AcfpParseError;
    fn try_from(input: &str) -> Result<Self, Self::Error> {
        let db_path: Vec<DotBracketVec> = input
            .split_whitespace()
            .map(|s| DotBracketVec::try_from(s).unwrap())
            .collect();
        let length = db_path.len();
        assert_eq!(length, db_path.last().unwrap().len());

        let pt_path: Vec<PairTable> = db_path
            .iter()
            .map(|db| PairTable::try_from(db).unwrap())
            .collect();

        let all_pairs: AHashSet<Pair> = pt_path
            .iter()
            .flat_map(|pt| {
                PairSet::from(pt).iter().collect::<Vec<Pair>>()
            })
            .collect();

        let mut po = PartialOrder::new();
        let mut partial_order = Some(po.clone());
        for npt in &pt_path {
            if !po.extend_by_pairtable(npt) {
                partial_order = None;
                break;
            }
        }
        if partial_order.is_some() {
            partial_order = Some(po);
        }

        Ok(Acfp { 
            db_path, 
            pt_path,
            all_pairs,
            partial_order,
            length,
        })
    }
}

pub fn enum_structs(n: usize) -> (usize, Vec<usize>) {
    let mut seq = vec![1, 1];

    for c in 2..=n {
        let prev = seq[c - 1];
        let mut next = 0;
        for i in 1..c {
            next += seq[c - i - 1] * seq[i - 1];
        }
        seq.push(prev + next);
    }
    (seq[n], seq)
}

pub fn generate_structs(n: usize) -> Vec<Vec<DotBracket>> {
    use DotBracket::*;
    fn helper(n: usize, memo: &mut HashMap<usize, Vec<Vec<DotBracket>>>
    ) -> Vec<Vec<DotBracket>> {
        if let Some(res) = memo.get(&n) {
            return res.clone();
        }

        let mut results = Vec::new();

        if n == 0 {
            results.push(vec![]);
        } else {
            // Add Unpaired to all structures of length n - 1
            for mut s in helper(n - 1, memo) {
                let mut new_s = vec![Unpaired];
                new_s.append(&mut s);
                results.push(new_s);
            }

            // Add all valid pairings: Open + inner + Close + outer
            for i in 1..n {
                let left = helper(i - 1, memo);
                let right = helper(n - i - 1, memo);
                for l in &left {
                    for r in &right {
                        let mut structure = vec![Open];
                        structure.extend_from_slice(l);
                        structure.push(Close);
                        structure.extend_from_slice(r);
                        results.push(structure);
                    }
                }
            }
        }

        memo.insert(n, results.clone());
        results
    }

    helper(n, &mut HashMap::new())
}

/// --- Core filtered cartesian product function ---
pub fn filtered_acfps<'a>(
    all_struct_sets: &'a [Vec<Vec<DotBracket>>],
) -> impl Iterator<Item = Acfp> + 'a {
    fn recurse<'a>(
        sets: &'a [Vec<Vec<DotBracket>>],    // remaining structure sets
        acfp: Acfp,
    ) -> Box<dyn Iterator<Item = Acfp> + 'a> {
        if !acfp.is_valid(&mut DomainRegistry::new()) {
            return Box::new(std::iter::empty());
        }

        if sets.is_empty() {
            return Box::new(std::iter::once(acfp));
        }

        let next_set = &sets[0];
        let rest = &sets[1..];

        Box::new(next_set.iter().flat_map(move |next| {
            let mut new_acfp = acfp.clone();
            new_acfp.extend_by_one(DotBracketVec(next.clone()));
            recurse(rest, new_acfp)
        }))
    }

    // Start recursion with an empty prefix
    recurse(all_struct_sets, Acfp::try_from(".").expect("must work"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acfp_01() {
        let acfp = Acfp::try_from(". () .() (()) ").unwrap();
        assert!(acfp.is_cycle_free());
        assert!(acfp.has_valid_pair_hierarchy());
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_acfp_02() {
        let acfp = Acfp::try_from(". () (.) (()) ").unwrap();
        assert!(!acfp.is_cycle_free());
        assert!(acfp.has_valid_pair_hierarchy());
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));

        let ph = acfp.pair_hierarchy().unwrap();
        assert_eq!(ph.get(&(0,1)), Some(&1));
        assert_eq!(ph.get(&(0,2)), Some(&2));
        assert_eq!(ph.get(&(0,3)), Some(&3));
        assert_eq!(ph.get(&(1,2)), Some(&1));
    }

    #[test]
    fn test_acfp_03() {
        let acfp = Acfp::try_from(". () (). ()() (()). ").unwrap();
        assert!(!acfp.is_cycle_free());
        assert!(!acfp.has_valid_pair_hierarchy());
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_acfp_04() {
        let acfp = Acfp::try_from(". () (). ()() (().) ").unwrap();
        println!("{:?}", acfp.pair_hierarchy());
        assert!(acfp.is_cycle_free());
        assert!(acfp.has_valid_pair_hierarchy());
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_acfp_05() {
        let acfp = Acfp::try_from(". .. ... .(.) ((.)) .(().)").unwrap();
        assert!(acfp.is_cycle_free());
        assert!(!acfp.has_valid_pair_hierarchy());
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_acfp_06() {
        let acfp = Acfp::try_from(". () (.) (.). (.)() ((..))").unwrap();
        assert!(acfp.is_cycle_free());
        assert!(acfp.has_valid_pair_hierarchy());
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_many_invalid_acfps() {
        let acfp = Acfp::try_from(". () (). ()() (.()) ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). ()() .(()) ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) ()() (.()) (.)(.)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) ()() .(()) (.)(.)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). (()) (()). ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). ().. (().) ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). .(.) ((.)) ().(.)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) (()). .((.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) (()). ()(..)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) ()(.) .((.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) ()(.) ()(.).").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) ()(.) ()(..)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) ()(.) (())()").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) ()(.) (.(.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() .(). (().) ()(..)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. .() (()) (()). .(().)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. .() (()) (()). .((.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) (.). (.)() ((.).)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() (()) (()). ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() .(). (().) ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_many_valid_acfps() {
        let acfp = Acfp::try_from(". () (.) (.). (.)() ((..))").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. (.) (..) (..). ((..))").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. (.) (()) (()). (()())").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) (.). ()(.) ()(())").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() .(). ()(.) ()(())").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_special_acfps() {
        let acfp = Acfp::try_from(". () .() .(). .()() .((.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) (.). (.)() (.(.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. (.) ..() (.)() ((.).)").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_implicit_invalid_acfps() {
        let acfp = Acfp::try_from(". .. .() ..() (.()) ((()))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. (.) (.). (().) (()())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). ().. (().) (()())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). (()) (()). (()())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (.) ()() ()(). ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () .() ()() ()(). ()(())").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_saturation_issues_acfps() {
        let acfp = Acfp::try_from(". .. .() (()) (()). (.(.))").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. .() (()) (()). ((..))").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". .. (.) (.). (.)() ((..))").unwrap();
        assert!(acfp.is_valid(&mut DomainRegistry::new()));
    }

    
    #[test]
    fn test_unsaturated_acfps() {
        let acfp = Acfp::try_from(". () (). .(.) .(..) .((.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));

        let acfp = Acfp::try_from(". () (). (..) (...) (.(.))").unwrap();
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }

    #[test]
    fn test_pseudoknotted_acfps() {
        let acfp = Acfp::try_from(". () (). .(.) .(..) .(...) .(...). .((...))").unwrap();
        assert!(acfp.has_valid_design(&mut DomainRegistry::new()));
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
        let acfp = Acfp::try_from(". () (). (..) (...) (....) (....). (.(...))").unwrap();
        assert!(acfp.has_valid_design(&mut DomainRegistry::new()));
        assert!(!acfp.is_valid(&mut DomainRegistry::new()));
    }


}


