//! Fuzz test: the free energy of a multi-stranded complex, with or without 
//! pseudoknots, must not depend on the order in which its strands are written.

use std::collections::BTreeSet;
use std::panic::{catch_unwind, AssertUnwindSafe};
use ff_structure::{MultiPairTable, MultiStruct};
use ff_energy::{Fuzzyfold, EnergyModel, Gap, NucleotideVec, Ring, RingDecomposition};

/// A tiny deterministic generator, to keep the test free of dependencies.
struct Lcg(u64);
impl Lcg {
    fn below(&mut self, m: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % m
    }
}

/// A complex: strands of bases, and pairs between (strand, index) positions.
struct Complex {
    strands: Vec<Vec<char>>,
    pairs: Vec<((usize, usize), (usize, usize))>,
}

impl Complex {
    fn random(rng: &mut Lcg) -> Option<Self> {
        let n = 2 + rng.below(3);
        let lens: Vec<usize> = (0..n).map(|_| 6 + rng.below(11)).collect();
        let mut strands: Vec<Vec<char>> = lens.iter()
            .map(|&l| (0..l).map(|_| ['A', 'C', 'G', 'U'][rng.below(4)]).collect())
            .collect();
        let mut used: Vec<Vec<bool>> = lens.iter().map(|&l| vec![false; l]).collect();
        let mut pairs = Vec::new();
        for _ in 0..(2 + rng.below(5)) {
            let (s, t) = (rng.below(n), rng.below(n));
            let h = 1 + rng.below(4);
            let (i, j) = (rng.below(lens[s]), rng.below(lens[t]));
            // Helix (s, i..i+h) with (t, j-h+1..=j), antiparallel.
            if i + h > lens[s] || j + 1 < h { continue; }
            // On one strand, at least 3 unpaired bases between the two halves.
            let (a, b) = (i + h - 1, j + 1 - h);
            if s == t && b < a + 4 { continue; }
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
        // Connected: all strands linked by pairs.
        let mut comp: Vec<usize> = (0..n).collect();
        let find = |c: &Vec<usize>, mut x: usize| { while c[x] != x { x = c[x]; } x };
        for &((s, _), (t, _)) in &pairs { let (a, b) = (find(&comp, s), find(&comp, t)); comp[a] = b; }
        let root = find(&comp, 0);
        (0..n).all(|s| find(&comp, s) == root).then_some(Complex { strands, pairs })
    }

    /// The sequence and pair table with strands in the given order.
    fn drawing(&self, order: &[usize]) -> (String, MultiPairTable) {
        let mut offset = vec![0; self.strands.len()];
        let mut seq = String::new();
        for (r, &s) in order.iter().enumerate() {
            if r > 0 { seq.push('+'); }
            offset[s] = seq.len();
            seq.extend(self.strands[s].iter());
        }
        let blank: String = seq.chars().map(|c| if c == '+' { '+' } else { '.' }).collect();
        let mut mpt = MultiPairTable::try_from(blank.as_str()).unwrap();
        for &((s, i), (t, j)) in &self.pairs {
            let (a, b) = (offset[s] + i, offset[t] + j);
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
fn strand_order_invariance() {
    let model = Fuzzyfold::default();
    let mut rng = Lcg(2026);
    let (mut tested, mut with_pk, mut failed) = (0, 0, Vec::new());

    while tested < 3000 {
        let Some(complex) = Complex::random(&mut rng) else { continue };
        let drawings: Vec<(String, MultiPairTable)> = permutations(complex.strands.len())
            .iter().map(|o| complex.drawing(o)).collect();
        tested += 1;
        let pk = drawings.iter().any(|(_, mpt)| mpt.rings().iter()
            .any(|r| matches!(r, Ring::Pseudoknot(_))));
        with_pk += pk as usize;

        let energies: Vec<Result<i32, String>> = drawings.iter().map(|(seq, mpt)| {
            let seq = NucleotideVec::try_from_rna(seq).unwrap();
            catch_unwind(AssertUnwindSafe(|| model.energy_of_structure(&seq, mpt)))
                .map_err(|_| "panic".to_string())
                .and_then(|r| r.map_err(|e| format!("{e:?}")))
        }).collect();
        if energies.iter().any(|e| e != &energies[0]) {
            failed.push((drawings, energies));
        }
    }
    println!("tested {tested} complexes ({with_pk} with pseudoknots in some order), {} failed", 
        failed.len());
    for (drawings, energies) in failed.iter().take(3) {
        for ((seq, mpt), e) in drawings.iter().zip(energies) {
            let db = ff_structure::DotBracketVec::try_from(mpt).unwrap();
            println!("  {seq}\n  {db} {e:?}");
        }
        println!();
    }
    assert!(failed.is_empty());
}

/// A base in strand coordinates: (strand, position in the strand).
type Base = (usize, usize);
/// A helix end as seen from a loop: the base the loop arrives at, and its partner.
type End = (Base, Base);
/// A flank of a helix end: the 5' flank (before the base the loop arrives 
/// at) or the 3' flank (after its partner).
type Flank = (End, bool);

/// The decomposition of a drawing, in strand coordinates and with the 
/// pseudoknots expanded into what they are made of. Everything that 
/// depends on the drawing is gone: indices, where the strands are cut, and 
/// whether a loop is the exterior or a joint exterior loop.
#[derive(Debug, Default, PartialEq)]
struct Physical {
    /// Closed loops (outside a pseudoknot, or interior loops of a strained stem):
    /// scaled or not, and their helix ends.
    loops: Vec<(bool, Vec<End>)>,
    /// The flanks of helix ends that are exterior-type.
    exterior: BTreeSet<Flank>,
    /// The flanks of helix ends and the unpaired bases in a pseudoloop.
    pseudoloop: BTreeSet<Flank>,
    unpaired: BTreeSet<Base>,
    /// The pairs of strained stems, and the number of groups of cycles.
    strained: BTreeSet<(Base, Base)>,
    knots: usize,
}

impl Physical {
    fn new(seq: &str, mpt: &MultiPairTable, order: &[usize]) -> Self {
        // Strand coordinates of each index (None for a strand break).
        let mut base = Vec::new();
        let (mut strand, mut pos) = (0, 0);
        for c in seq.chars() {
            if c == '+' {
                base.push(None);
                strand += 1;
                pos = 0;
            } else {
                base.push(Some((order[strand], pos)));
                pos += 1;
            }
        }
        let at = |x: usize| base[x].unwrap();
        // A helix end as seen from a loop: (the base it arrives at, its partner).
        let end = |(a, b): (usize, usize)| (at(a), at(b));
        let pair = |(a, b): (usize, usize)| (at(a).min(at(b)), at(a).max(at(b)));
        // Both flanks of a helix end in an exterior-type loop.
        let both = |p: (usize, usize)| [(end(p), false), (end(p), true)];

        let mut phys = Physical::default();
        for r in mpt.rings() {
            match &r {
                Ring::Loop(l) => match l.gap() {
                    Gap::None => phys.loops.push((false, l.pairs().iter().map(|&p| end(p)).collect())),
                    Gap::Open { .. } => phys.exterior.extend(l.pairs().iter().flat_map(|&p| both(p))),
                    Gap::Wrap { .. } => unreachable!("No circles here."),
                },
                Ring::Pseudoknot(pk) => {
                    for l in pk.stacks.iter().chain(&pk.interior_loops) {
                        phys.loops.push((true, l.pairs().iter().map(|&p| end(p)).collect()));
                        phys.strained.extend(l.pairs().iter().map(|&p| pair(p)));
                    }
                    phys.knots += pk.knots;
                    for l in &pk.loops {
                        let pairs = l.nested.pairs();
                        let start = match l.nested.gap() {
                            Gap::Open { five, .. } => five,
                            _ => pairs[pairs.len() - 1].1 + 1,
                        };
                        let k = l.strained.len();
                        for (t, &(a, b)) in pairs.iter().enumerate() {
                            // The 5' flank lies in the segment before, the 3' flank in the one after.
                            for (three, strained) in [(false, l.strained[t]), (true, l.strained[(t + 1) % k])] {
                                if strained { &mut phys.pseudoloop } else { &mut phys.exterior }
                                    .insert((end((a, b)), three));
                            }
                            if l.stem_ends[t] {
                                phys.strained.insert(pair((a, b)));
                            }
                            if l.strained[t] {
                                let c = if t == 0 { start } else { pairs[t - 1].1 + 1 };
                                phys.unpaired.extend((c..a).map(at));
                            }
                        }
                    }
                }
            }
        }
        for (_, ends) in &mut phys.loops {
            ends.sort_unstable();
        }
        phys.loops.sort_unstable();
        phys
    }
}

#[test]
fn strand_order_invariance_of_the_decomposition() {
    let mut rng = Lcg(2026);
    let (mut tested, mut failed) = (0, Vec::new());
    while tested < 3000 {
        let Some(complex) = Complex::random(&mut rng) else { continue };
        tested += 1;
        let orders = permutations(complex.strands.len());
        let physical: Vec<(String, MultiPairTable, Physical)> = orders.iter().map(|o| {
            let (seq, mpt) = complex.drawing(o);
            let phys = Physical::new(&seq, &mpt, o);
            (seq, mpt, phys)
        }).collect();
        if physical.iter().any(|(_, _, p)| p != &physical[0].2) {
            failed.push(physical);
        }
    }
    println!("tested {tested} complexes, {} decompose differently", failed.len());
    for physical in failed.iter().take(2) {
        for (seq, mpt, phys) in physical {
            let db = ff_structure::DotBracketVec::try_from(mpt).unwrap();
            println!("  {seq}\n  {db}\n  {phys:?}");
        }
        println!();
    }
    assert!(failed.is_empty());
}

impl Complex {
    /// Two complexes side by side, without pairs between them.
    fn beside(&self, other: &Complex) -> Complex {
        let n = self.strands.len();
        let mut strands = self.strands.clone();
        strands.extend(other.strands.iter().cloned());
        let mut pairs = self.pairs.clone();
        pairs.extend(other.pairs.iter().map(|&((s, i), (t, j))| ((s + n, i), (t + n, j))));
        Complex { strands, pairs }
    }

    fn energy(&self, model: &Fuzzyfold, order: &[usize]) -> i32 {
        let (seq, mpt) = self.drawing(order);
        model.energy_of_structure(&NucleotideVec::try_from_rna(&seq).unwrap(), &mpt).unwrap()
    }
}

/// A random order of n strands.
fn shuffle(rng: &mut Lcg, n: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    for k in (1..n).rev() {
        order.swap(k, rng.below(k + 1));
    }
    order
}

#[test]
fn detached_complexes() {
    // Two complexes (or a complex and a single strand) drawn together, in 
    // random strand orders, have the energy of the two complexes apart.
    let model = Fuzzyfold::default();
    let mut rng = Lcg(2027);
    let (mut tested, mut failed) = (0, Vec::new());
    while tested < 1000 {
        let Some(a) = Complex::random(&mut rng) else { continue };
        let b = if rng.below(4) == 0 {
            let len = 4 + rng.below(8);
            Complex { strands: vec![(0..len).map(|_| ['A', 'C', 'G', 'U'][rng.below(4)]).collect()], pairs: vec![] }
        } else {
            let Some(b) = Complex::random(&mut rng) else { continue };
            b
        };
        tested += 1;
        let apart = a.energy(&model, &(0..a.strands.len()).collect::<Vec<_>>()) 
            + b.energy(&model, &(0..b.strands.len()).collect::<Vec<_>>());
        let both = a.beside(&b);
        for _ in 0..6 {
            let order = shuffle(&mut rng, both.strands.len());
            let together = catch_unwind(AssertUnwindSafe(|| both.energy(&model, &order)));
            if together.as_ref().ok() != Some(&apart) {
                let (seq, mpt) = both.drawing(&order);
                let db = ff_structure::DotBracketVec::try_from(&mpt).unwrap();
                failed.push(format!("{seq}\n  {db} {together:?}, apart {apart}"));
            }
        }
    }
    println!("tested {tested} pairs of complexes, {} drawings failed", failed.len());
    for f in failed.iter().take(5) {
        println!("  {f}");
    }
    assert!(failed.is_empty());
}
