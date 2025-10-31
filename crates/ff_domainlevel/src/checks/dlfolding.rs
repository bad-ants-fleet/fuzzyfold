use ff_structure::NAIDX;
use ndarray::Array2;
use rustc_hash::FxHashSet;
use rustc_hash::FxHashMap;
use crate::domain::{DomainRefVec, DomainRegistry};

pub fn nussinov(p: &Array2<usize>) -> Array2<usize> {
    let (n, m) = p.dim();
    assert!(n == m);
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
pub fn build_pair_scores(
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

#[allow(dead_code)] // for MFE prediction, no subopts.
fn traceback(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
    pairs: &mut Vec<(usize, usize)>,
) {
    if i >= j {
        return;
    }
    let dp_ij = dp[(i,j)];

    if dp_ij == dp[(i + 1, j)] {
        traceback(i + 1, j, dp, p, pairs);
    } else if dp_ij == dp[(i, j - 1)] {
        traceback(i, j - 1, dp, p, pairs);
    } else if p[(i, j)] > 0 && dp_ij == dp[(i + 1, j - 1)] + p[(i, j)] {
        pairs.push((i, j));
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

pub fn traceback_structures(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
) -> Vec<Vec<Option<NAIDX>>> {
 
    let mut memo = FxHashMap::default();
    let as_pairs = traceback_all(i, j, dp, p, &mut memo);

    as_pairs.into_iter()
        .map(|ps| {
            let mut table = vec![None; j + 1];
            for (p, q) in ps {
                table[p] = Some(q as NAIDX);
                table[q] = Some(p as NAIDX);
            }
            table
        })
        .collect()
}

pub fn traceback_all(
    i: usize,
    j: usize,
    dp: &Array2<usize>,
    p: &Array2<usize>,
    memo: &mut FxHashMap<(usize, usize), FxHashSet<Vec<(usize, usize)>>>,
) -> FxHashSet<Vec<(usize, usize)>> {
    if i >= j {
        return FxHashSet::from([vec![]].into_iter().collect());
    }

    if let Some(cached) = memo.get(&(i, j)) {
        return cached.clone();
    }

    let mut results = FxHashSet::default();
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
            assert!(i < j);
            sub.push((i, j));
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
    memo.insert((i, j), results.clone());
    results
}


#[cfg(test)]
mod tests {
    use super::*;

    fn drv(input_seq: &str, reg: &DomainRegistry) -> DomainRefVec {
        input_seq
            .split_whitespace()
            .map(|name| reg.get(name).unwrap())
            .collect()
    }

    #[test]
    fn test_pair_score_simple() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("b", 1);
        registry.intern("c", 1);

        let input = "a a* b b* c";
        let domains = drv(input, &registry);
        let p = build_pair_scores(&domains, &registry);
        assert_eq!(p[(0, 1)], 1);
        assert_eq!(p[(1, 0)], 1);
        assert_eq!(p[(2, 3)], 1);
        assert_eq!(p[(3, 2)], 1);
        assert_eq!(p[(0, 2)], 0);
    }

    #[test]
    fn test_nussinov_basic_structure() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("b", 2);
 
        let domains = drv("a a* b b*", &registry);
        let p = build_pair_scores(&domains, &registry);
        let dp = nussinov(&p);
        assert_eq!(dp[(0, 3)], 3); // a-a* and b-b*

        let mut pairs: Vec<(usize, usize)> = Vec::new();
        traceback(0, domains.len() - 1, &dp, &p, &mut pairs);
        assert_eq!(pairs, vec![(0, 1), (2, 3)]); // a-a* and b-b*
    }

    #[test]
    fn test_traceback_all_variants() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        registry.intern("x", 2);
        let sequence = drv("a x a*", &registry);
        let p = build_pair_scores(&sequence, &registry);
        assert_eq!(p[(0, 2)], 1); // ensure complement
        assert_eq!(p[(2, 0)], 1); // ensure complement
        let dp = nussinov(&p);
        let mut memo = FxHashMap::default();
        let structs: Vec<Vec<(usize, usize)>> = traceback_all(0, sequence.len() - 1, &dp, &p, &mut memo).into_iter().collect();
        assert_eq!(dp[(0, 2)], 1);
        assert_eq!(structs.len(), 1);
        assert!(structs[0].contains(&(0, 2)));
    }

    #[test]
    fn test_traceback_all_bifurcation() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        let sequence = drv("a a* a a*", &registry);
        let p = build_pair_scores(&sequence, &registry);
        let dp = nussinov(&p);
        let mut memo = FxHashMap::default();
        let all_structs: Vec<Vec<(usize, usize)>> = traceback_all(0, sequence.len() - 1, &dp, &p, &mut memo).into_iter().collect();
        println!("{:?}", all_structs);
        assert_eq!(dp[(0, 3)], 2);

        assert!(all_structs.iter().any(|s| s.contains(&(0, 1)) && s.contains(&(2, 3))));
        assert!(all_structs.iter().any(|s| s.contains(&(0, 3)) && s.contains(&(1, 2))));
        let mut all_structs = traceback_structures(0, sequence.len() - 1, &dp, &p);
        all_structs.sort_unstable();
        assert_eq!(all_structs, [[Some(1), Some(0), Some(3), Some(2)],
                                 [Some(3), Some(2), Some(1), Some(0)]]);
    }

    #[test]
    fn test_traceback_all_multioutput() {
        let mut registry = DomainRegistry::new();
        registry.intern("a", 1);
        let sequence = drv("a a* a a* a a* a a*", &registry);
        let p = build_pair_scores(&sequence, &registry);
        let dp = nussinov(&p);
        let all_structs: Vec<Vec<(usize, usize)>> = traceback_all(0, sequence.len() - 1, &dp, &p, &mut FxHashMap::default()).into_iter().collect();
        println!("{:?}", all_structs);
        assert_eq!(all_structs.len(), 14);
        let all_structs = traceback_structures(0, sequence.len() - 1, &dp, &p);
        for s in &all_structs {
            println!("{:?}", s);
        }
        assert_eq!(all_structs.len(), 14);
    }
}

