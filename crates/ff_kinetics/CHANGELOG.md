# Changelog ff_kinetics

All notable changes to this crate will be documented in this file.

## [development]
### Fixed
- SSA::simulate applied the move whose time fell after t_max: a simulation
  ended in the structure after it, and in co-transcriptional simulations
  one extra move happened before each extension. Output recorded at fixed
  times was not affected. The simulation now ends in the structure of the
  last callback (also the Python SimulationIterator). Simulations that end
  at t_max draw one random number less, so seeded runs after the first
  differ (the no-shifts sum of kinetics_workload is now -5305625).
- The neighborhood enumeration (enum_neighbors::enumerate_neighbors, was
  LoopNeighbors::generate_neighbors): with a step limit, a structure first
  found in more moves than needed was not expanded again, and structures
  within the limit were missed; with shift moves, the energy budget was
  reduced by the barriers along the path instead of the energy. The budget
  is now the energy of each structure, and a test compares the result with
  a breadth-first search. It works on any MoveSet and returns it to the
  start. ff-explore and the Python Explorer round the energy delta.
- Co-transcriptional simulations with shift moves: after an extension, the
  shift moves of the exterior loop kept their old rates and missed the new
  base. A test compares the moves after every step and extension with a
  move set built from scratch.

### Added
- Crate documentation: move sets, rate models, the SSA, co-transcriptional
  folding, output times, macrostates and time courses, neighborhoods.
- OutputTimes: the times at which a simulation reports its state; any
  increasing list, lin_log (the scheme of the command line tools),
  transcription (one per extension, then logarithmic), and a Default (50
  times up to 1e-4 s, 50 in log space up to 1 s).
- SSA::sample: simulate (co-transcriptionally) and record the move set at
  each output time; extension_schedule gives the times of co_simulate.
  ff-timecourse, ff-ensemble and the Python bindings use them.

### Changed
- K0 (0 degrees Celsius in Kelvin) is called ZERO_CELSIUS, not to be
  confused with the rate constant k0.
- Timecourse::plot_svg returns a Result instead of panicking.
- MacrostateRegistry::insert_files reads each file once (not once per
  length of a co-transcriptional macrostate), and ignores anything after
  a structure (such as its energy) without a warning; the Python bindings
  use it instead of their own reader.
- Tests assert instead of printing; docstrings ASCII only.
- Building a move set (LoopTable, LoopNeighbors) fails with MoveSetError
  (Unsupported, Energy) instead of a String.
- Only the move set changes its loop table: the mutating and lookup
  functions of LoopTable, and the apply_*_move functions and neighbor maps
  of LoopNeighbors, are crate-internal. The modules loop_table,
  three_way_shifts and four_way_shifts are private (LoopTable is exported
  from the crate root).
- rate_tree is private to the SSA, documented, and built in one step
  (FromIterator instead of init_insert and init_partial_sums); dump_by_rate
  and the unused alternative update of the SSA are removed. Tests compare
  the tree with a map under random operations, and the moves and rates
  after every step and extension with a move set built from scratch.
  Identical trajectories.
- The modules timeline, timeline_io and timeline_plotting are one module,
  timecourse: `Timecourse { times, counts }` holds the counts of each
  macrostate (registry index) at each output time, without the registry,
  and is not generic over the energy model. The caller classifies
  structures (MacrostateRegistry::classify) and adds a trajectory, one
  macrostate per output time (add_trajectory); merge adds the counts.
  The registry is an argument where names are needed: read_tln,
  write_tln, write_nxy (was Display) and plot_svg (was
  plot_occupancy_over_time). The .tln and .nxy formats are unchanged.
- criterion is a dev-dependency; ndarray is no longer a dependency.
- A structure in two macrostates of the same length is an error when the
  macrostate files are read (MacrostateRegistry::insert_from_reader),
  instead of a panic in MacrostateRegistry::classify.
- The trait Walker is called MoveSet, and its lists of moves with their
  energy changes (Moves) are called Neighbors.
- MoveSet is only what a simulation needs (propose and apply moves, the
  current structure and energy); the trait Elongation adds the extension
  at the 3' end, and SSA::co_simulate requires it.
- SSA::waiting_time and SSA::apply_random_move are the two halves of a
  step, SSA::apply_extension an extension; simulate and co_simulate use
  them, with identical trajectories. SSA::move_set gives the move set.
- Shift moves are a runtime setting, `ShiftMoves { three_way, four_way }`,
  instead of the type parameter `ShiftPolicy` (NoShift, ThreeWayOnly,
  FourWayOnly, ThreeAndFour): `LoopNeighbors<E>` takes it on construction,
  and `Arrhenius::shift_moves()` gives the moves that have a rate. Identical
  trajectories, within 0.3% of the instructions.
- LoopNeighbors: the moves after a change come from two functions, the
  delete moves of the pairs of the new loops (with the shift moves of the
  loops on their other side) and the moves of the new loops, instead of
  six copies. Identical trajectories, within 0.7% of the instructions.
  Walker documents the moves that apply_move returns.
- Indices are usize (Move, LoopTable, LoopNeighbors), not NAIDX
- LoopTable reads the loops of the decomposition and rejects structures
  with crossing pairs
- LoopTable::try_from returns an error, not a panic, for a structure that
  cannot be evaluated (e.g. a hairpin that is too small)

### Added
- The move sets work on NestedLoop, the loops of the ring decomposition
  (NearestNeighborLoop is removed): identical trajectories, about 4% (no
  shifts) and 8% (shifts) more instructions. Crossing pairs are a
  pseudoknot also where they close no cycle
- benches/workloads/kinetics_workload.rs (an example target): seeded
  simulations for instruction counts
- tests/trajectories.rs: digests of seeded trajectories with and without
  shifts, on linear and circular strands; at every step the flux is the
  sum of the rates of a move set built from scratch
- Circular strands, alone and at full length: add, delete and shift moves
  through the strand wrap; moves are proposed only if their loops fit
  (LoopTable::fits: the hairpin through the wrap of a circle is not too
  small). Structures are written with `^`

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

