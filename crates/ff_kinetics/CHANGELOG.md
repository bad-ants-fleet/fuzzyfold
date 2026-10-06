# Changelog ff_kinetics

All notable changes to this crate will be documented in this file.

## 0.5.0 - 2026-10-06
### Changed
- The same macrostate can now be defined at multiple sequence lengths. A flag
allows for automatic generation of macrostates for shorter lengths during
co-transcriptional folding. 
- Timeline printing and visualization has been adapted to merge macrostates 
with the same name across different sequence lengths.
- co_simulate does now modify the callback to reduce tinc in case it goes beyond
the current transcription step. 

### Removed
- preliminary interface to merge timelines.

## 0.4.2 - 2026-06-18
### Changed
- no more lifetimes for macrostate interface to facilitate exports to python.
- default timeline visualization.

### Added
- preliminary interface to merge timelines.
- shift_analysis feature to report shift-move activation energies as they are calculated.

## 0.4.1 - 2026-03-04
### Changed
- added Arrhenius, removed Metropolis/Kawasaki
- removed lifetimes and added Arc<EnergyModel> for Python exports.
- LoopTable now owns the sequence! 
- compatibility-updates to ff_energy v0.4

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

