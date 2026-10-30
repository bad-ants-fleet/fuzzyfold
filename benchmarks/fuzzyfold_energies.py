#!/usr/bin/env python3
"""Recorded energies of the pseudoknotted benchmark sets of ff_energy.

For each structure of benchmark_random_pseudoknots_len{50,100}.vrna, the
energy of fuzzyfold's default model (Fuzzyfold, rna_extended, 37 C) in its
four variants: one strand, four strands, circular, and four strands with the
third circular (cut at the quarters), in units of 10 cal/mol. The regression
test tests/pseudoknot_energies.rs compares with them; record them again
after an intended change of the model.

Needs the fuzzyfold Python package (crates/ff_python).

    python benchmarks/fuzzyfold_energies.py
"""

import sys
from pathlib import Path

import fuzzyfold as ff

DATA = Path(__file__).resolve().parent.parent / "crates" / "ff_energy" / "benches" / "data"
LENGTHS = [50, 100]
VARIANTS = ["one", "four", "circular", "four-circular"]


def quarters(s, ends):
    """The string cut into four at the quarters, with the given strand ends."""
    n = len(s)
    cuts = [0, n // 4, n // 2, 3 * n // 4, n]
    parts = [s[cuts[k]:cuts[k + 1]] for k in range(4)]
    return parts[0] + ends[0] + parts[1] + ends[1] + parts[2] + ends[2] + parts[3]


def variant(s, name):
    return {
        "one": s,
        "four": quarters(s, "+++"),
        "circular": s + "^",
        "four-circular": quarters(s, "++^"),
    }[name]


def main():
    model = ff.Fuzzyfold()
    for length in LENGTHS:
        lines = (DATA / f"benchmark_random_pseudoknots_len{length}.vrna").read_text().splitlines()
        out = ["# Fuzzyfold (rna_extended, 37 C) energies in 10 cal/mol of",
               f"# benchmark_random_pseudoknots_len{length}.vrna, one line per structure:",
               "# " + " ".join(VARIANTS)]
        for k in range(0, len(lines), 3):
            seq, db = lines[k + 1], lines[k + 2]
            out.append(" ".join(
                str(round(model.energy_of_structure(variant(seq, v), variant(db, v)) * 100))
                for v in VARIANTS))
        path = DATA / f"benchmark_random_pseudoknots_len{length}_fuzzyfold.txt"
        path.write_text("\n".join(out) + "\n")
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
