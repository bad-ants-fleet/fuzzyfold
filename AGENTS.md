# Notes for coding agents

fuzzyfold is a Rust workspace for nucleic acid secondary structures: energy
evaluation (with pseudoknots and multiple strands) and stochastic folding
kinetics. One author reviews every line; read this before changing code.

## Crates

| Crate          | Owns                                                         | Depends on               |
|----------------|--------------------------------------------------------------|--------------------------|
| `ff_structure` | structure representations, closed regions, bands, density    | -                        |
| `ff_energy`    | sequences, loop decomposition, nearest neighbor energy models | `ff_structure`           |
| `ff_kinetics`  | move sets, rate models, stochastic simulation, macrostates   | `ff_structure`, `ff_energy` |
| `fuzzyfold`    | command line tools (`src/bin/ff-*.rs`), re-exports the above  | all                      |
| `ff_python`    | Python bindings (PyO3, built with maturin)                   | all                      |

## How the pieces fit

- Structures: `DotBracketVec` (1 byte per position, for storage and hashing)
  and `MultiPairTable` (2 bytes per position, the input of everything else).
  Both hold one or more strands and pseudoknots. See the crate docs of
  `ff_structure` (`cargo doc -p ff_structure --open`).
- Energy: `MultiPairTable` implements `RingDecomposition`, which yields
  `Ring` items: a `NestedLoop` (a loop traced along the backbone) or a
  `Pseudoknot`; an `EnergyModel` evaluates them: `Fuzzyfold` (under
  development) or `ViennaRNA` (fixed, for comparisons; `Fuzzyfold` uses its
  loop energies until the models diverge). A `Pseudoknot` is not a loop: it
  holds the loops on cycles of the loop graph, a block (biconnected
  component) of loops joined by pairs, and which segments are strained
  (`decomposition/pseudoknots.rs`); `Fuzzyfold::pseudoknot()`
  (`models/fuzzyfold.rs`) evaluates it with `PKParams`
  (`parameters/pk_parameterset.rs`).
- Kinetics: `LoopTable` is built once from a `MultiPairTable`, with the
  `NestedLoop`s of its ring decomposition; `LoopNeighbors` proposes moves
  and updates the table; `SSA` simulates. `tests/trajectories.rs` pins
  seeded trajectories (with flux and energy checks at every step).
  Kinetics support one strand without pseudoknots, linear or circular (at
  full length), checked in `LoopTable::try_from`.

## Conventions

- Indices are positions in the structure string; the end of a strand takes
  one position: a strand break (`+` or `&`) or a strand wrap (`^`, the
  strand before it is circular). Indices are `usize`; only
  `MultiPairTable` stores positions in two bytes (at most 65533).
- Energies are `i32` in units of 10 cal/mol (100 = 1 kcal/mol); the Python
  bindings return kcal/mol.
- Bracket levels: `()` is 0, `[]` 1, `{}` 2, `<>` 3, `Aa`..`Zz` 4..29.
  Writing a structure assigns levels again, from `(` upward.

## Working rules

- Discuss the design before writing code, and make no code changes without
  explicit permission.
- Simple, explainable changes. No premature optimization: do not add small
  changes of uncertain benefit.
- Docstrings short, ASCII only, no decorative symbols. No attribution or
  author names in code.
- Prefer established terms from the literature over new ones.
- Keep the CHANGELOG of each changed crate up to date (`## [development]`,
  followed by the version of the next release, which is already set in
  `Cargo.toml`). All crates share one version; the first breaking change
  after a release raises the minor version (0.x).

## Checking work

```bash
cargo test --workspace                 # unit tests, doctests, fuzz tests
cargo clippy --workspace
cargo bench -p ff_energy -p ff_kinetics  # about 5 minutes
```

Timings on a shared machine vary by 10% or more between runs. To compare
performance, count instructions instead on the fixed, seeded workloads, and
check that the results they print (sums of energies) are identical:

```bash
cargo build --release -p ff_energy -p ff_kinetics --examples
valgrind --tool=callgrind target/release/examples/energy_workload pseudoknots four
valgrind --tool=callgrind target/release/examples/kinetics_workload shifts
```

The workloads are in `benches/workloads` of each crate (built as examples,
for a fixed path). The benchmark sets are in `crates/ff_energy/benches/data`:
nested structures (`structures`) and pseudoknotted ones (`pseudoknots`),
each evaluated as one strand, four strands, circular, and four strands with
the third circular. The scripts that make them, and the runtime
comparisons with other programs, are in `benchmarks/`.

Python bindings: `cd crates/ff_python && maturin develop --release`, then
`pytest tests`.

## Before merging into `development`

Go through this list for every branch, also when the changes were made
with an LLM. Each item names what to check and how.

Structure
- Code is in the crate that owns it (table above). Library crates have no
  command line dependencies (clap, anyhow, colored, indicatif); argument
  parsing lives in `fuzzyfold`, whose binaries import the libraries
  through `fuzzyfold::structure`, `fuzzyfold::energy` and
  `fuzzyfold::kinetics`.
- Layout of a crate: unit tests next to the code; `tests/` for integration
  tests (recorded data in `tests/data`); `benches/` for criterion benches
  and `benches/workloads/*_workload.rs` for instruction counts. Scripts
  that need other programs (ViennaRNA, Python packages) are in
  `benchmarks/`.
- Dependencies are declared in `[workspace.dependencies]` and used by name
  (`rand.workspace = true`). No unused dependencies: grep each one in the
  crate that lists it. `cargo search <name>` for current versions.
- No stray files: `*.drf`, `*.tln`, `*.nxy`, `*.so` and `__pycache__` are
  not committed. `Cargo.lock` is committed (its diffs are hidden by
  `.gitattributes`); it changes with the dependencies, not on its own.

Correctness
- `cargo test --workspace --all-features` passes (the feature
  `shift_analysis` compiles code the default build does not).
- `tests/trajectories.rs` pins seeded trajectories. Record them again
  (`FF_RECORD_TRAJECTORIES=1`) only for an intended change of the order of
  moves or random numbers, and say why in the CHANGELOG. To show that only
  the order changed, compare both versions with the moves of
  `propose_moves` sorted: the workload sums must be identical.
- A change to a move set keeps the test that compares the rate tree with a
  move set built from scratch after every move and extension.

Performance
- Count instructions (callgrind, see above) before and after, on the
  workloads with the same results. Report the numbers in the commit
  message or the CHANGELOG. Keep a change for speed only if it is
  measurable and simple; drop it otherwise.
- On hot paths, prefer vectors indexed by loop index or position over hash
  maps, and plain iteration over collect-and-clone. The rate tree index is
  the one hash map there.

Readability
- Simple, explainable code; no wrappers or helpers that hide a design
  problem instead of fixing it. Established terms; names that say what a
  thing is.
- Docstrings short, ASCII only, no attribution. Comments explain choices
  that are not obvious, not what the next line does.
- `cargo clippy --workspace --all-targets --all-features` reports nothing.

Documentation
- Public entry points have a short, realistic example: no seeds, a comment
  on every value that is not obvious (rate constants and their units,
  times, energy units, what a callback returns), output printed rather
  than asserted unless it is exact (energies).
- `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps` builds and
  all doctests pass.
- The `## [development]` section of the CHANGELOG of each changed crate
  says what changed for a user, including changed seeded outputs and
  workload sums.
- `README.md`, the READMEs of the crates (their examples run as doctests),
  the READMEs in `examples/` and `TODO.md` match the code (command line
  options, file names, open items).
- This file still describes the crates, how they fit, and how to check
  work.

Interfaces (command line, Python)
- A change to the Rust API is carried through to the command line tools
  (`kinetics_parsers.rs`, `energy_parsers.rs`) and the Python bindings
  (`crates/ff_python`), with their help texts and CHANGELOGs.
- Python: `maturin develop --release`, then `pytest tests` in
  `crates/ff_python` (energies, structures, trajectories, ensembles, time
  courses, edits, neighborhoods) and `examples/example.py`. The bindings use the stable
  ABI (`abi3-py310`): one wheel per platform for CPython 3.10 and later.

Process
- Design discussed before code; one topic per commit, with the reason in
  the message. Work stays on its `dev_*` branch until this list is done.

## Branches

`main` is stable, `development` integrates features, `dev_*` branches hold
work in progress. Open items are in `TODO.md`.
