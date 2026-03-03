# The fuzzyfold workspace
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![fuzzyfold](https://img.shields.io/crates/v/fuzzyfold.svg?label=fuzzyfold)](https://crates.io/crates/fuzzyfold)
[![ff_structure](https://img.shields.io/crates/v/ff_structure.svg?label=ff_structure)](https://crates.io/crates/ff_structure)
[![ff_energy](https://img.shields.io/crates/v/ff_energy.svg?label=ff_energy)](https://crates.io/crates/ff_energy)
[![ff_kinetics](https://img.shields.io/crates/v/ff_kinetics.svg?label=ff_kinetics)](https://crates.io/crates/ff_kinetics)
[![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)]()

An open-source collection of nucleic acid folding algorithms.

**Note**: This is a _very_ early stage, rapidly developing coding project. You
are welcome to use it for research, but be prepared for frustration from
drastic interface changes. You may use GitHub issues for suggestions, but you
are also welcome to reach out directly at this point.

## Current fuzzyfold software
 - **ff-eval**: Free-energy evaluation for secondary structures.
 - **ff-trajectory**: Single stochastic nucleic acid folding trajectories.
 - **ff-timecourse**: Stochastic nucleic acid secondary structure ensemble simulations.
 - **ff-randseq**: Generate a random sequence. 
 - **ff-explore**: Enumerate secondary structure neighborhoods.

(Other software is work in progress and not yet published to crates.io)

## Current fuzzyfold crates
 - ff_structure: Nucleic acid secondary structure data structures.
 - ff_energy: Secondary structure free energy evaluation.
 - ff_kinetics: Stochastic folding kinetics for nucleic acids.

(Other crates are work in progress and not yet published to crates.io)

## Developer notes
Thank you for considering contributing! The goal of the fuzzyfold workspace is
to provide a coherent, well-documented ecosystem for RNA and DNA
secondary-structure modeling, kinetic simulations, and analysis. Each crate
focuses on a clearly separated aspect of the workflow: structure, energy
evaluation, and/or kinetic modeling. 

We welcome improvements of any kind, from bug fixes and performance
enhancements to documentation, examples, and new features. Feel free to reach
out with specific ideas.

For benchmarking of the stochastic simulation algorithm:

```cargo bench --workspace```

### Git branches (work in progress)
 - **main**: well-structured, well-documented, high-coverage, publication-ready code.
 - **development**: well-integrated, some documentation, some coverage, experimental code.
 - **dev_feature**: early code/crate proposals.



