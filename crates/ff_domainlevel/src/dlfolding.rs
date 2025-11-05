//! Nussinov-style domain-level base-pair maximization.
//!

use ndarray::Array2;
use ahash::AHashSet;
use nohash_hasher::IntMap;

use ff_structure::NAIDX;
use ff_structure::P1KEY;
use ff_structure::Pair;
use ff_structure::PairSet;
use ff_structure::DotBracket;
use ff_structure::DotBracketVec;

use crate::DomainRefVec;
use crate::DomainRegistry;
use crate::error::RegistryError;


/// Nussinov-style domain-level folding algorithm.
/// Can be initialized from a domain sequence and the corresponding registry to
/// look up the domain lengths, or through an arbitrary matrix of pair scores.
pub struct NussinovDP {
    pair_scores: Array2<usize>,
    dp_table: Array2<usize>,
}

impl From<Array2<usize>> for NussinovDP {
    fn from(pair_scores: Array2<usize>) -> Self {
        let dp_table = nussinov(&pair_scores);
        Self {
            pair_scores,
            dp_table,
        }
    }
}

impl From<(&DomainRefVec, &DomainRegistry)> for NussinovDP {
    fn from((domains, registry): (&DomainRefVec, &DomainRegistry)) -> Self {
        let pair_scores = build_pair_scores(domains, registry);
        let dp_table = nussinov(&pair_scores);
        Self {
            pair_scores,
            dp_table,
        }
    }
}

impl TryFrom<(&str, &DomainRegistry)> for NussinovDP {
    type Error = RegistryError;

    fn try_from((sequence, registry): (&str, &DomainRegistry)) -> Result<Self, Self::Error> {
        let domains: Result<DomainRefVec, _> = sequence
            .split_whitespace()
            .map(|name| registry.get(name)
                .ok_or(RegistryError::UnknownDomain(name.to_string())))
            .collect();

        let domains = domains?;
        Ok(NussinovDP::from((&domains, registry)))
    }
}

impl NussinovDP {

    pub fn get_mfe_pairs(&self, len: Option<usize>) -> PairSet {
        let length = match len {
            Some(l) => l,
            None => self.dp_table.nrows(),
        };
        let mut pairs = PairSet::new(length);
        traceback(0, length - 1, &self.dp_table, &self.pair_scores, &mut pairs);
        pairs
    }

    pub fn all_mfe_pairs(&self, len: Option<usize>) -> Vec<PairSet> {
        let length = match len {
            Some(l) => l,
            None => self.dp_table.nrows(),
        };
        let mut memo: IntMap<P1KEY, AHashSet<Vec<P1KEY>>> = IntMap::default();
        let as_pairs = traceback_all(0, length - 1, &self.dp_table, &self.pair_scores, &mut memo);
        as_pairs.into_iter()
            .map(|ps| {
                let mut pset = PairSet::new(length);
                for p1key in ps {
                    pset.insert(Pair::from_key(p1key));
                }
                pset
            })
        .collect()
    }

    pub fn all_mfe_structs(&self, len: Option<usize>) -> Vec<DotBracketVec> {
        let length = match len {
            Some(l) => l,
            None => self.dp_table.nrows(),
        };
        let mut memo: IntMap<P1KEY, AHashSet<Vec<P1KEY>>> = IntMap::default();
        let as_pairs = traceback_all(0, length - 1, &self.dp_table, &self.pair_scores, &mut memo);
        as_pairs.into_iter()
            .map(|ps| {
                let mut dbv = vec![DotBracket::Unpaired; length];
                for p1key in ps {
                    let pair = Pair::from_key(p1key);
                    dbv[pair.i() as usize] = DotBracket::Open;
                    dbv[pair.j() as usize] = DotBracket::Close;
                }
                DotBracketVec(dbv)
            })
        .collect()
    }

    pub fn pair_scores(&self) -> &Array2<usize> {
        &self.pair_scores
    }

    pub fn dp_table(&self) -> &Array2<usize> {
        &self.dp_table
    }
}


fn nussinov(p: &Array2<usize>) -> Array2<usize> {
    let (n, m) = p.dim();
    debug_assert!(n == m);
    let mut dp = Array2::from_elem((n, n), 0);
    for l in 1..n {
        for i in 0..n - l {
            let j = i + l;
            let mut max_val = dp[(i + 1, j)].max(dp[(i, j - 1)]);
            if p[(i, j)] > 0 {
                max_val = max_val.max(dp[(i + 1, j - 1)] + p[(i, j)]);
            }
            for k in i + 1..j {
                max_val = max_val.max(dp[(i, k)] + dp[(k + 1, j)]);
            }
            dp[(i, j)] = max_val;
        }
    }
    dp
}

/// Returns a pairwise score matrix for a vector of Domains.
fn build_pair_scores(
    domains: &DomainRefVec, 
    registry: &DomainRegistry
) -> Array2<usize> {
    let n = domains.len();
    let mut p = Array2::from_elem((n, n), 0);

    for ((i, j), value) in p.indexed_iter_mut() {
        assert_eq!(*value, 0); // sanity check
        let di = &domains[i];
        let dj = &domains[j];
        if registry.are_complements(di, dj) {
            *value = di.length.min(dj.length);
        }
    }
    p
}

fn traceback(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
    pairs: &mut PairSet,
) {
    if i >= j {
        return;
    }
    let dp_ij = dp[(i, j)];

    if dp_ij == dp[(i + 1, j)] {
        traceback(i + 1, j, dp, p, pairs);
    } else if dp_ij == dp[(i, j - 1)] {
        traceback(i, j - 1, dp, p, pairs);
    } else if p[(i, j)] > 0 && dp_ij == dp[(i + 1, j - 1)] + p[(i, j)] {
        pairs.insert(Pair::new(i as NAIDX, j as NAIDX));
        traceback(i + 1, j - 1, dp, p, pairs);
    } else {
        for k in i + 1..j {
            if dp_ij == dp[(i, k)] + dp[(k + 1, j)] {
                traceback(i, k, dp, p, pairs);
                traceback(k + 1, j, dp, p, pairs);
                break;
            }
        }
    }
}

fn traceback_all(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
    memo: &mut IntMap<P1KEY, AHashSet<Vec<P1KEY>>>,
) -> AHashSet<Vec<P1KEY>> {
    if i >= j {
        return AHashSet::from([vec![]]);
    }

    if let Some(cached) = memo.get(&(Pair::new(i as NAIDX, j as NAIDX).key())) {
        return cached.clone();
    }

    let mut results = AHashSet::default();
    let dp_ij = dp[(i,j)];

    // Case 1: i unpaired
    if dp_ij == dp[(i + 1, j)] {
        for sub in traceback_all(i + 1, j, dp, p, memo) {
            results.insert(sub);
        }
    // Case 2: j unpaired 
    } else if dp_ij == dp[(i, j - 1)] {
        for sub in traceback_all(i, j - 1, dp, p, memo) {
            results.insert(sub);
        }
    }

    // Case 3: i-j paired
    if p[(i, j)] > 0 && dp_ij == dp[(i + 1, j - 1)] + p[(i, j)] {
        for mut sub in traceback_all(i + 1, j - 1, dp, p, memo) {
            debug_assert!(i < j);
            sub.push(Pair::new(i as NAIDX, j as NAIDX).key());
            sub.sort_unstable();
            results.insert(sub);
        }
    }

    // Case 4: bifurcation
    for k in i + 1..j {
        if dp_ij == dp[(i, k)] + dp[(k + 1, j)] {
            let lefts = traceback_all(i, k, dp, p, memo);
            let rights = traceback_all(k + 1, j, dp, p, memo);

            if lefts.is_empty() || rights.is_empty() {
                continue;
            }

            for left in &lefts {
                for right in &rights {
                    let mut combined = left.clone();
                    combined.extend(right);
                    results.insert(combined);
                }
            }
        }
    }
    memo.insert(Pair::new(i as NAIDX, j as NAIDX).key(), results.clone());
    results
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_score_simple() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("b", 1);
        registry.intern("c", 1);

        let ndp = NussinovDP::try_from(("a a* b b* c", &registry)).unwrap();
        assert_eq!(ndp.pair_scores[(0, 1)], 1);
        assert_eq!(ndp.pair_scores[(1, 0)], 1);
        assert_eq!(ndp.pair_scores[(2, 3)], 1);
        assert_eq!(ndp.pair_scores[(3, 2)], 1);
        assert_eq!(ndp.pair_scores[(0, 2)], 0);
    }

    #[test]
    fn test_nussinov_basic_structure() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("b", 2);

        let ndp = NussinovDP::try_from(("a a* b b*", &registry)).unwrap();
        assert_eq!(ndp.dp_table[(0, 3)], 3);
        let pairs = ndp.get_mfe_pairs(None);
        assert!(pairs.contains(&Pair::new(0, 1)));
        assert!(pairs.contains(&Pair::new(2, 3)));
    }

    #[test]
    fn test_traceback_all_variants() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("x", 2);

        let ndp = NussinovDP::try_from(("a x a*", &registry)).unwrap();
        assert_eq!(ndp.pair_scores[(0, 2)], 1);
        assert_eq!(ndp.pair_scores[(2, 0)], 1);
        assert_eq!(ndp.dp_table[(0, 2)], 1);

        let structs = ndp.all_mfe_structs(None);
        assert_eq!(structs.len(), 1);
        assert_eq!(structs[0], DotBracketVec::try_from("(.)").unwrap());
    }

    #[test]
    fn test_traceback_all_bifurcation() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);

        let ndp = NussinovDP::try_from(("a a* a a*", &registry)).unwrap();
        assert_eq!(ndp.dp_table[(0, 3)], 2);

        let structs = ndp.all_mfe_structs(None);
        assert_eq!(structs.len(), 2);
        assert!(structs.contains(&DotBracketVec::try_from("(())").unwrap()));
        assert!(structs.contains(&DotBracketVec::try_from("()()").unwrap()));
    }

    #[test]
    fn test_traceback_all_multioutput() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        let ndp = NussinovDP::try_from(("a a* a a* a a* a a*", &registry)).unwrap();
        let structs = ndp.all_mfe_structs(None);
        assert_eq!(structs.len(), 14);
    }
}

