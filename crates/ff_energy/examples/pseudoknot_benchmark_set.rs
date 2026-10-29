//! Generate the pseudoknotted benchmark sets from the nested ones: each
//! structure is cut into segments that are written in another order, with
//! every pair kept between its (moved) bases, so that pairs between the
//! segments cross. A third of the structures each get a pseudoknot density
//! (Rastegari & Condon, `ClosedRegions::density`) of 2, 3, and 4 or more:
//! cuts at random into 2, 3, ... segments and shuffles them (seeded) until
//! the density fits. Pairs that could close a hairpin with fewer than three
//! unpaired bases are dropped: inside, or through the wrap of the circular
//! variants (the whole strand, and its third quarter). Every variant of the
//! benchmark (one strand, four strands, circular, four strands with a
//! circular third) must evaluate.
//!
//! cargo run --release -p ff_energy --example pseudoknot_benchmark_set

use std::fmt::Write;

use ff_structure::{ClosedRegions, DotBracketVec, MultiPairTable, MultiStruct};
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/benches/data");
const LENGTHS: [usize; 7] = [50, 100, 250, 500, 750, 1000, 2500];

/// The string cut into four at the quarters, with the given strand ends.
fn quarters(s: &str, ends: [char; 3]) -> String {
    let n = s.len();
    let cuts = [n / 4, n / 2, 3 * n / 4];
    let mut out = String::with_capacity(n + 3);
    let mut from = 0;
    for (cut, end) in cuts.into_iter().zip(ends) {
        out.push_str(&s[from..cut]);
        out.push(end);
        from = cut;
    }
    out.push_str(&s[from..]);
    out
}

/// The variants of a structure in the benchmark: one strand, four strands,
/// circular, and four strands with the third circular.
fn variants(seq: &str, db: &str) -> [(String, String); 4] {
    [
        (seq.to_string(), db.to_string()),
        (quarters(seq, ['+', '+', '+']), quarters(db, ['+', '+', '+'])),
        (format!("{seq}^"), format!("{db}^")),
        (quarters(seq, ['+', '+', '^']), quarters(db, ['+', '+', '^'])),
    ]
}

/// A tiny deterministic generator, to keep the example free of dependencies.
struct Lcg(u64);
impl Lcg {
    fn below(&mut self, m: usize) -> usize {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % m
    }
}

/// The structure cut at `cuts` (ascending, from 0 to its length) and
/// written in the given order of segments, with every pair kept between
/// its (moved) bases, except those that could close a hairpin of fewer than
/// three bases (inside, or through the wrap of a circular variant).
fn permuted(seq: &str, db: &str, cuts: &[usize], order: &[usize]) -> (String, MultiPairTable) {
    let n = seq.len();
    let start = |s: usize| cuts[s];
    // The new position of each base.
    let mut pos = vec![0; n];
    let mut new_seq = String::with_capacity(n);
    for &s in order {
        for (i, p) in pos.iter_mut().enumerate().take(start(s + 1)).skip(start(s)) {
            *p = new_seq.len() + (i - start(s));
        }
        new_seq.push_str(&seq[start(s)..start(s + 1)]);
    }
    // The circles of the variants, (first, last): the strand, its third quarter.
    let circles = [(0, n - 1), (n / 2, 3 * n / 4 - 1)];
    let fits = |p: usize, q: usize| q - p > 3 && circles.iter()
        .all(|&(f, l)| !(f <= p && q <= l) || (l - q) + (p - f) >= 3);
    let old = MultiPairTable::try_from(db).unwrap();
    let mut new = MultiPairTable::try_from(".".repeat(n).as_str()).unwrap();
    for i in 0..n {
        if let MultiStruct::Paired(j) = old.get(i) && i < j {
            let (p, q) = (pos[i].min(pos[j]), pos[i].max(pos[j]));
            if fits(p, q) {
                new.set(p, MultiStruct::Paired(q));
                new.set(q, MultiStruct::Paired(p));
            }
        }
    }
    (new_seq, new)
}

/// Whether a density belongs to a class: 2, 3, or 4 and more.
fn in_class(density: usize, class: usize) -> bool {
    if class < 4 { density == class } else { density >= 4 }
}

fn main() {
    let model = Fuzzyfold::default();
    for len in LENGTHS {
        let input = std::fs::read_to_string(format!("{DATA}/benchmark_random_structures_len{len}.vrna")).unwrap();
        let lines: Vec<&str> = input.lines().collect();
        let mut out = String::new();
        let mut found = [0; 3];
        let mut densities = std::collections::BTreeMap::new();
        for (r, record) in lines.chunks(3).enumerate() {
            let [name, seq, db] = record else { panic!("Malformed input.") };
            // A third of the structures each of density 2, 3, and 4 or more:
            // cut into more and more segments at random and shuffle them,
            // until the density fits.
            let class = 2 + r % 3;
            let mut rng = Lcg(r as u64);
            let mut best = None;
            'search: for k in 2..=12 {
                for _ in 0..50 {
                    let mut order: Vec<usize> = (0..k).collect();
                    for t in (1..k).rev() {
                        order.swap(t, rng.below(t + 1));
                    }
                    // Cut points at random, every segment at least one base.
                    let mut cuts: Vec<usize> = (0..k - 1).map(|_| 1 + rng.below(seq.len() - 1)).collect();
                    cuts.extend([0, seq.len()]);
                    cuts.sort_unstable();
                    cuts.dedup();
                    if cuts.len() != k + 1 {
                        continue;
                    }
                    let (s, mpt) = permuted(seq, db, &cuts, &order);
                    let density = ClosedRegions::from(&mpt).density();
                    if in_class(density, class) {
                        best = Some((s, mpt));
                        found[class - 2] += 1;
                        *densities.entry(density).or_insert(0) += 1;
                        break 'search;
                    }
                }
            }
            let Some((seq, mpt)) = best else {
                println!("{name}: no structure of density class {class}");
                continue;
            };
            let db = DotBracketVec::try_from(&mpt).expect("Too many pseudoknot levels.").to_string();
            let evaluates = variants(&seq, &db).iter().all(|(s, d)| {
                let (Ok(s), Ok(d)) = (NucleotideVec::try_from(s.as_str()), MultiPairTable::try_from(d.as_str())) else {
                    return false;
                };
                model.energy_of_structure(&s, &d).is_ok()
            });
            assert!(evaluates, "{name}: a variant does not evaluate.");
            writeln!(out, "{name}_pk\n{seq}\n{db}").unwrap();
        }
        std::fs::write(format!("{DATA}/benchmark_random_pseudoknots_len{len}.vrna"), out).unwrap();
        println!("length {len}: density 2, 3, 4+: {found:?}, densities {densities:?}");
    }
}
