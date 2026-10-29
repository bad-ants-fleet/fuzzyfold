//! Fuzz test: circular strands together with linear strands, or with each
//! other. The energy must not depend on the order of the strands, nor on
//! where the circles are cut.

use ff_structure::{MultiPairTable, MultiStruct};
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};

/// A tiny deterministic generator, to keep the test free of dependencies.
struct Lcg(u64);
impl Lcg {
    fn below(&mut self, m: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % m
    }
}

/// Strands, some circular, and pairs between (strand, index) positions.
struct Complex {
    strands: Vec<Vec<char>>,
    circular: Vec<bool>,
    pairs: Vec<((usize, usize), (usize, usize))>,
}

impl Complex {
    fn random(rng: &mut Lcg) -> Self {
        let n = 2 + rng.below(2);
        // One or two circles, the rest linear.
        let circles = 1 + rng.below(2);
        let circular: Vec<bool> = (0..n).map(|s| s < circles).collect();
        let lens: Vec<usize> = (0..n).map(|_| 8 + rng.below(12)).collect();
        let mut strands: Vec<Vec<char>> = lens.iter()
            .map(|&l| (0..l).map(|_| ['A', 'C', 'G', 'U'][rng.below(4)]).collect())
            .collect();
        let mut used: Vec<Vec<bool>> = lens.iter().map(|&l| vec![false; l]).collect();
        let mut pairs = Vec::new();
        for _ in 0..(1 + rng.below(5)) {
            let (s, t) = (rng.below(n), rng.below(n));
            let h = 1 + rng.below(4);
            let (i, j) = (rng.below(lens[s]), rng.below(lens[t]));
            if i + h > lens[s] || j + 1 < h { continue; }
            // On one strand: at least 3 unpaired bases between the innermost pair.
            if s == t && j + 1 - h < i + h + 3 { continue; }
            if (0..h).any(|d| used[s][i + d] || used[t][j - d]) { continue; }
            for d in 0..h {
                used[s][i + d] = true;
                used[t][j - d] = true;
                let (a, b) = [('G', 'C'), ('C', 'G'), ('A', 'U'), ('U', 'A')][rng.below(4)];
                strands[s][i + d] = a;
                strands[t][j - d] = b;
                pairs.push(((s, i + d), (t, j - d)));
            }
        }
        Complex { strands, circular, pairs }
    }

    /// The sequence and pair table with the strands in `order`, each circle
    /// `s` cut `cut[s]` bases later.
    fn drawing(&self, order: &[usize], cut: &[usize]) -> (String, MultiPairTable) {
        let len = |s: usize| self.strands[s].len();
        let local = |s: usize, i: usize| if self.circular[s] { (i + len(s) - cut[s]) % len(s) } else { i };
        let mut offset = vec![0; self.strands.len()];
        let mut seq = String::new();
        for (k, &s) in order.iter().enumerate() {
            if k > 0 && !seq.ends_with('^') { seq.push('+'); }
            offset[s] = seq.len();
            seq.extend((0..len(s)).map(|i| self.strands[s][if self.circular[s] { (i + cut[s]) % len(s) } else { i }]));
            if self.circular[s] { seq.push('^'); }
        }
        let blank: String = seq.chars().map(|c| if "+^".contains(c) { c } else { '.' }).collect();
        let mut mpt = MultiPairTable::try_from(blank.as_str()).unwrap();
        for &((s, i), (t, j)) in &self.pairs {
            let (a, b) = (offset[s] + local(s, i), offset[t] + local(t, j));
            mpt.set(a, MultiStruct::Paired(b));
            mpt.set(b, MultiStruct::Paired(a));
        }
        (seq, mpt)
    }
}

fn permutations(n: usize) -> Vec<Vec<usize>> {
    if n == 1 { return vec![vec![0]]; }
    let mut out = Vec::new();
    for p in permutations(n - 1) {
        for k in 0..n {
            let mut q = p.clone();
            q.insert(k, n - 1);
            out.push(q);
        }
    }
    out
}

#[test]
fn circles_with_linear_strands() {
    let model = Fuzzyfold::default();
    let mut rng = Lcg(2028);
    let (mut tested, mut joined) = (0, 0);
    while tested < 1000 {
        let complex = Complex::random(&mut rng);
        tested += 1;
        joined += complex.pairs.iter().any(|&((s, _), (t, _))| s != t) as usize;
        let n = complex.strands.len();
        let mut energies = Vec::new();
        for order in permutations(n) {
            for random in [false, true] {
                let cut: Vec<usize> = (0..n).map(|s| match random && complex.circular[s] {
                    true => rng.below(complex.strands[s].len()),
                    false => 0,
                }).collect();
                let (seq, mpt) = complex.drawing(&order, &cut);
                let e = model.energy_of_structure(&NucleotideVec::try_from_rna(&seq).unwrap(), &mpt)
                    .map_err(|e| format!("{e:?}"));
                energies.push((seq, ff_structure::DotBracketVec::try_from(&mpt).unwrap().to_string(), e));
            }
        }
        for (seq, db, e) in &energies {
            assert_eq!(e, &energies[0].2, "\\n{}\\n{} {:?}\\n{seq}\\n{db} {e:?}", energies[0].0, energies[0].1, energies[0].2);
        }
    }
    assert!(joined > tested / 2, "{joined} of {tested} with strands joined");
}
