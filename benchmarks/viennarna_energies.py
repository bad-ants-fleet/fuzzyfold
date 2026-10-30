#!/usr/bin/env python3
"""Reference energies of the nested benchmark sets of ff_energy from ViennaRNA.

For each structure of benchmark_random_structures_len{50,100}.vrna, the
energy in ViennaRNA (eval_structure, dangles 2, as RNAeval -d2, at 37 C)
with the parameter sets Turner 2004, Andronescu 2007 and Mathews 2004 (DNA;
ViennaRNA treats U as T), in units of 10 cal/mol. The model ViennaRNA of
ff_energy must reproduce them (tests/reference_energies.rs).

Needs the ViennaRNA Python bindings (pip install ViennaRNA).

    python benchmarks/viennarna_energies.py
"""

import sys
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "crates" / "ff_energy" / "benches" / "data"
LENGTHS = [50, 100]
PARAMETERS = ["turner2004", "andronescu2007", "mathews2004"]


def records(length):
    lines = (DATA / f"benchmark_random_structures_len{length}.vrna").read_text().splitlines()
    return [(lines[k + 1], lines[k + 2]) for k in range(0, len(lines), 3)]


def energies(params):
    """The energies of all sets with one parameter set. Run in a process of
    its own: loading a parameter set on top of another may keep parts of it."""
    import RNA
    load = {
        "turner2004": RNA.params_load_RNA_Turner2004,
        "andronescu2007": RNA.params_load_RNA_Andronescu2007,
        "mathews2004": RNA.params_load_DNA_Mathews2004,
    }
    load[params]()
    md = RNA.md()
    md.dangles = 2
    md.temperature = 37.0
    return {length: [round(RNA.fold_compound(seq, md).eval_structure(db) * 100)
                     for seq, db in records(length)]
            for length in LENGTHS}


def main():
    import RNA
    results = {}
    for params in PARAMETERS:
        with ProcessPoolExecutor(1, mp_context=get_context("spawn")) as pool:
            results[params] = pool.submit(energies, params).result()
    for length in LENGTHS:
        out = [f"# ViennaRNA {RNA.__version__} energies (as RNAeval -d2, 37 C) in 10 cal/mol of",
               f"# benchmark_random_structures_len{length}.vrna, one line per structure:",
               "# " + " ".join(PARAMETERS)]
        for k in range(len(records(length))):
            out.append(" ".join(str(results[p][length][k]) for p in PARAMETERS))
        path = DATA / f"benchmark_random_structures_len{length}_viennarna.txt"
        path.write_text("\n".join(out) + "\n")
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
