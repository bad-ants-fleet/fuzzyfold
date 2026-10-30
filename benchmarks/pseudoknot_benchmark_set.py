#!/usr/bin/env python3
"""Generate the pseudoknotted benchmark sets of ff_energy from the nested ones.

Each structure is cut into segments that are written in another order, with
every pair kept between its (moved) bases, so that pairs between the segments
cross. A third of the structures each get a pseudoknot density (Rastegari &
Condon, fuzzyfold.density) of 2, 3, and 4 or more: cut at random into 2, 3,
... segments and shuffled (seeded) until the density fits. Pairs that could
close a hairpin with fewer than three unpaired bases are dropped: inside, or
through the wrap of the circular variants (the whole strand, and its third
quarter). Every variant of the benchmark (one strand, four strands,
circular, four strands with a circular third) must evaluate.

Needs the fuzzyfold Python package (crates/ff_python).

    python benchmarks/pseudoknot_benchmark_set.py
"""

import sys
from collections import Counter
from pathlib import Path

import fuzzyfold as ff

DATA = Path(__file__).resolve().parent.parent / "crates" / "ff_energy" / "benches" / "data"
LENGTHS = [50, 100, 250, 500, 750, 1000, 2500]
BRACKETS = ["()", "[]", "{}", "<>"] + [c + c.lower() for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"]
MASK = (1 << 64) - 1


class Lcg:
    """A small deterministic generator (the one of the former Rust version)."""

    def __init__(self, seed):
        self.state = seed

    def below(self, m):
        self.state = (self.state * 6364136223846793005 + 1442695040888963407) & MASK
        return (self.state >> 33) % m


def quarters(s, ends):
    """The string cut into four at the quarters, with the given strand ends."""
    n = len(s)
    cuts = [0, n // 4, n // 2, 3 * n // 4, n]
    parts = [s[cuts[k]:cuts[k + 1]] for k in range(4)]
    return parts[0] + ends[0] + parts[1] + ends[1] + parts[2] + ends[2] + parts[3]


def variants(seq, db):
    """One strand, four strands, circular, four strands with the third circular."""
    return [
        (seq, db),
        (quarters(seq, "+++"), quarters(db, "+++")),
        (seq + "^", db + "^"),
        (quarters(seq, "++^"), quarters(db, "++^")),
    ]


def dot_bracket(n, pairs):
    """A dot-bracket string of the pairs, normalized (levels from '(' up)."""
    levels = []
    chars = ["."] * n
    for i, j in sorted(pairs):
        for level, taken in enumerate(levels):
            if all(not (a < i < b < j or i < a < j < b) for a, b in taken):
                taken.append((i, j))
                break
        else:
            level = len(levels)
            levels.append([(i, j)])
        chars[i], chars[j] = BRACKETS[level][0], BRACKETS[level][1]
    return ff.normalize("".join(chars))


def permuted(seq, db, cuts, order):
    """The structure cut at cuts and written in the given order of segments,
    with every pair kept between its (moved) bases, except those that could
    close a hairpin of fewer than three bases (inside, or through the wrap
    of a circular variant)."""
    n = len(seq)
    pos = [0] * n
    new_seq = ""
    for s in order:
        for i in range(cuts[s], cuts[s + 1]):
            pos[i] = len(new_seq) + i - cuts[s]
        new_seq += seq[cuts[s]:cuts[s + 1]]
    circles = [(0, n - 1), (n // 2, 3 * n // 4 - 1)]

    def fits(p, q):
        return q - p > 3 and all(not (f <= p and q <= l) or (l - q) + (p - f) >= 3
                                 for f, l in circles)

    pairs = []
    for i, j in ff.pairs(db):
        p, q = min(pos[i], pos[j]), max(pos[i], pos[j])
        if fits(p, q):
            pairs.append((p, q))
    return new_seq, pairs


def in_class(density, cls):
    """Whether a density belongs to a class: 2, 3, or 4 and more."""
    return density == cls if cls < 4 else density >= 4


def main():
    model = ff.Fuzzyfold()
    for length in LENGTHS:
        lines = (DATA / f"benchmark_random_structures_len{length}.vrna").read_text().splitlines()
        out = []
        found = [0, 0, 0]
        densities = Counter()
        for r in range(len(lines) // 3):
            name, seq, db = lines[3 * r:3 * r + 3]
            cls = 2 + r % 3
            rng = Lcg(r)
            best = None
            for k in range(2, 13):
                for _ in range(50):
                    order = list(range(k))
                    for t in range(k - 1, 0, -1):
                        u = rng.below(t + 1)
                        order[t], order[u] = order[u], order[t]
                    # Cut points at random, every segment at least one base.
                    cuts = sorted({1 + rng.below(len(seq) - 1) for _ in range(k - 1)} | {0, len(seq)})
                    if len(cuts) != k + 1:
                        continue
                    new_seq, pairs = permuted(seq, db, cuts, order)
                    new_db = dot_bracket(len(new_seq), pairs)
                    density = ff.density(new_db)
                    if in_class(density, cls):
                        best = (new_seq, new_db)
                        found[cls - 2] += 1
                        densities[density] += 1
                        break
                if best:
                    break
            if best is None:
                print(f"{name}: no structure of density class {cls}")
                continue
            for s, d in variants(*best):
                model.energy_of_structure(s, d)  # raises if a variant does not evaluate
            out.append(f"{name}_pk\n{best[0]}\n{best[1]}\n")
        (DATA / f"benchmark_random_pseudoknots_len{length}.vrna").write_text("".join(out))
        print(f"length {length}: density 2, 3, 4+: {found}, densities {dict(sorted(densities.items()))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
