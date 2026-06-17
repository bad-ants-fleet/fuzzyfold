# fuzzyfold

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/bad-ants-fleet/fuzzyfold/blob/development/LICENSE)
[![fuzzyfold](https://img.shields.io/crates/v/fuzzyfold.svg?label=fuzzyfold)](https://crates.io/crates/fuzzyfold)
[![ff_structure](https://img.shields.io/crates/v/ff_structure.svg?label=ff_structure)](https://crates.io/crates/ff_structure)
[![ff_energy](https://img.shields.io/crates/v/ff_energy.svg?label=ff_energy)](https://crates.io/crates/ff_energy)
[![ff_kinetics](https://img.shields.io/crates/v/ff_kinetics.svg?label=ff_kinetics)](https://crates.io/crates/ff_kinetics)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](https://github.com/bad-ants-fleet/fuzzyfold/blob/development)

A high-performance framework for stochastic RNA folding kinetics.

`fuzzyfold` is an open-source Rust workspace for nucleic acid secondary
structure analysis with an explicit focus on kinetic modeling. It combines
a Gillespie-type stochastic simulation engine with a flexible move set 
(including three-way and four-way shift moves) and provides both command-line
tools for end users and a library interface with Python bindings for custom
workflows.

## Features

- Stochastic simulation of RNA and DNA secondary structure folding at
  base-pair resolution
- Three-way and four-way shift moves with Arrhenius-based activation energy
  parameterization
- Ensemble-level occupancy analysis via aggregation of parallel trajectories
- Co-transcriptional folding trajectories
- Thermodynamic parameters: RNA (Turner 2004), DNA (RNAstructure), extended
  sets including modification parameters and special hairpin parameters
- Energy evaluation fully consistent with ViennaRNA
- Pure Rust workspace: straightforward installation via `cargo`, reproducible
  builds, no external runtime dependencies
- Python bindings for rapid prototyping and custom analysis workflows

## Installation

`fuzzyfold` requires a working [Rust toolchain](https://rustup.rs). To install
all command-line tools:

```bash
cargo install fuzzyfold
```

Individual crates can also be installed separately:

```bash
cargo install ff_structure
cargo install ff_energy
cargo install ff_kinetics
```

## Command-line tools

| Program | Description |
|---|---|
| `ff-trajectory` | Single stochastic folding trajectory |
| `ff-timecourse` | Ensemble occupancy analysis over multiple parallel trajectories |
| `ff-explore` | Enumerate secondary structure neighborhoods and macrostates |
| `ff-eval` | Free-energy evaluation for secondary structures |
| `ff-randseq` | Generate a random nucleic acid sequence |

See directory [examples] for basic usage of command-line software,
as well as [examples/workflows] to reproduce published analyses.

## Crates

| Crate | Description |
|---|---|
| `ff_structure` | Nucleic acid secondary structure data structures |
| `ff_energy` | Free-energy evaluation using nearest-neighbor parameters |
| `ff_kinetics` | Stochastic folding kinetics engine |
| `fuzzyfold` | Command-line interfaces |

Additional crates are in development and will be published to
[crates.io](https://crates.io) as they mature.

## Python bindings

Python bindings are available for all major simulation and analysis components,
enabling rapid prototyping before contributing high-performance Rust code. See
the [crates/ff_python] directory for examples and installation instructions.

## Repository structure

| Branch | Description |
|---|---|
| `main` | Stable, well-documented, publication-ready code |
| `development` | Integrated development branch; some experimental features |
| `dev_feature` | Early crate or feature proposals |

## Contributing

Contributions are welcome! Please start with small contributions, such as bug
fixes, performance improvements, documentation, or examples. For new features
or new crates associated with fuzzyfold, please reach out before submitting a
pull request. The goal of the workspace is a coherent, well-documented
ecosystem for RNA and DNA secondary-structure modeling and kinetic simulation,
with each crate covering a clearly separated aspect of the workflow.

To run the benchmark suite:

```bash
cargo bench --workspace
```

Please open a GitHub issue to discuss larger changes before submitting a pull
request.

## Citation

If you use `fuzzyfold` in your research, please cite:

Stefan Badelt: **fuzzyfold: a high-performance framework for stochastic RNA folding kinetics**, 2026,
[https://doi.org/]

## License

MIT — see [LICENSE] for details.
