# Changelog ff_python

All notable changes to this crate will be documented in this file.

## [development] - 0.6.0
### Added
- Package metadata: readme, author, license, keywords, classifiers and
  URLs (repository, documentation, citation).
- The Python package fuzzyfold around the extension module
  (fuzzyfold._fuzzyfold), with numpy as a dependency.
- fuzzyfold.Timecourse: the macrostate counts of many simulations at the
  output times; occupancy, standard_error, merge, read (.tln) and save
  (.tln, .nxy, the formats of ff-timecourse).
- fuzzyfold.plotting (with matplotlib): plot_timecourse and the command
  line `python -m fuzzyfold.plotting` (was examples/py-utils/plot_anxy.py),
  with standard error bands; --t-split is called --t-sep and is given in
  the units of the input, like --t-ext.

- lin_log_times and transcription_times: output times for simulations.

- A-to-I editing: Simulator(k_atoi=...), SimulationIterator.sequence, the
  sequence in the keys of simulate_ensemble and as a column of .drf files
  (with edits), and macrostates of edited sequences.

### Changed
- simulate_ensemble writes the .drf file with ff_kinetics::ensemble: the
  time column is exact (it was rounded to three decimals).
- pyo3 0.30 and rand 0.10 (random generators seeded with
  rand::make_rng); serde_json is no longer a dependency.
- Tests with pytest (tests/test_fuzzyfold.py; the extra `test` installs
  pytest).
- Built against the stable ABI of CPython 3.10 (abi3-py310): one wheel
  for all CPython versions from 3.10 on, also those newer than pyo3.
  The extension-module feature is set in Cargo.toml only.
- Simulator.simulate_timecourse returns a Timecourse.
- simulate_ensemble and simulate_timecourse take a list of output times
  (`times`) instead of t_end, t_sep, t_lin and t_log; the simulation ends
  at the last one. The default is lin_log_times() (to 1 s), or with t_ext
  transcription_times (to 1 s after the full length).
- The bindings no longer depend on the fuzzyfold (command line) crate.

### Removed
- Plotter (use fuzzyfold.plotting) and
  Simulator.aggregate_timecourse_results (use Timecourse.merge, or
  `python -m fuzzyfold.plotting` with several .tln files and --merged).


## 0.5.0 - 2026-10-06
### Changed
- Default rna parameter specification.
- simulate -> simulate_trajectory.

### Added
- export for simulate_ensemble.
- export for simulate_timecourse.
- various file i/o conversions
