# Benchmark utilities

Scripts that make the benchmark sets of the crates and the reference energies
the tests compare with, and runtime comparisons of the command line tools with
other programs. Only these scripts need other programs (ViennaRNA, Python). The benchmarks themselves
(criterion, and the workloads for instruction counts) are in `benches/` of
each crate.

- `gen_benchmarkset.sh`: random sequences (`ff-randseq`) and their minimum
  free energy structures (`RNAsubopt -p 1`): the nested benchmark sets
  (`benchmark_random_structures_len*.vrna`, in `crates/ff_energy/benches/data`).
  Writes into the current directory.
- `pseudoknot_benchmark_set.py`: the pseudoknotted sets
  (`benchmark_random_pseudoknots_len*.vrna`) from the nested ones, in
  `crates/ff_energy/benches/data`. Needs the fuzzyfold Python package;
  deterministic (about 6 minutes).
- `viennarna_energies.py`: reference energies of the nested sets (lengths
  50 and 100) from ViennaRNA, for `tests/reference_energies.rs` of ff_energy.
  Needs the ViennaRNA Python bindings.
- `fuzzyfold_energies.py`: the energies of the pseudoknotted sets (lengths
  50 and 100, all variants) for the regression test
  `tests/pseudoknot_energies.rs`; record them again after an intended change
  of the model. Needs the fuzzyfold Python package.
- `eval_benchmark.sh`, `simulate_benchmark.sh`, `cotr_benchmark.sh`:
  wall-clock runtimes of `ff-eval`, `ff-trajectory` and other programs
  (`RNAeval`, `Kinfold`) on the benchmark sets, written as CSV;
  `plot_single_runtimes.py` plots them.
