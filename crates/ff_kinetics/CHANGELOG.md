# Changelog ff_kinetics

All notable changes to this crate will be documented in this file.

## [0.4.0] - 2026-02-06
### Added
- three-way and four-way shift moves for Metropolis model.
- co-transcriptional simulations.

### Changed
- removed LoopStructure, replaced it by:
   LoopTable, LoopNeighbors, and the Walker trait.
- renamed explore.rs to enum_neighbors.rs.
- new plotting style (lin/log equal split).

## [0.3.0] - 2026-01-13
### Added
- rate_tree.rs: logarithmic neighbor selection in stochastic simulations (huge speedup).
- explore.rs: basic neighborhood (macrostate) exploration (ff-locmin).
- rate_model.rs: Kawasaki rate model.

### Changed
- Major rewriting of LoopStructure and LoopStructureSSA
- Using FxHashMap instead of AHashMap for performance.
- Cargo bench updates.
- Macrostate naming.

### Removed
- Reaction for SSA.
- Preliminary implementations for commit & delay.

