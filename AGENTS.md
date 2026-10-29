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

`crates/ff_findpath` is not a workspace member and does not build against the
current crates. Ignore it unless asked.

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
- Keep the CHANGELOG of each changed crate up to date (`## [development]`).

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

The benchmark sets are in `crates/ff_energy/benches/data`: nested
structures (`structures`) and pseudoknotted ones (`pseudoknots`, made from
them by `examples/pseudoknot_benchmark_set.rs`), each evaluated as one
strand, four strands, circular, and four strands with the third circular.

Python bindings: `cd crates/ff_python && maturin develop --release`.

## Branches

`main` is stable, `development` integrates features, `dev_*` branches hold
work in progress. Open items are in `TODO.md`.
