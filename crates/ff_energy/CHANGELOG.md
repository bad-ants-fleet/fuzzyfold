# Changelog ff_energy

All notable changes to this crate will be documented in this file.

## [development] - 0.6.0
### Fixed
- A stack of a strained stem through the wraps of two circles was scaled
  as an interior loop: whether two pairs stack was decided by index. The
  decomposition now decides it along the backbone, and `Pseudoknot` holds
  the `stacks` of its strained stems apart from its `interior_loops`.
  With different scales for stacks and interior loops, the energy then
  depended on where the circles were cut. The energy terms (`ff-eval -v`)
  name these loops `Interior` with their pairs, also through wraps, where
  they were listed as `Pseudol.`.
- At temperatures other than 37 C, the extrapolation of loops longer than
  30 bases (lxc) was multiplied by the temperature in Celsius instead of
  the ratio of absolute temperatures: a hairpin of 40 bases at 25 C was
  7.5 kcal/mol too unfavorable. Tests against RNAeval (Turner 2004,
  Andronescu 2007, Mathews 2004)

### Added
- A README of the crate, with an example tested with the doctests.
- Benchmarks of pseudoknotted structures: benchmark_random_pseudoknots_len*
  (each nested structure cut into segments written in another order, a
  third each of pseudoknot density 2, 3, and 4 or more; without the pairs
  that could close a hairpin under three bases, also through a wrap; made
  by benchmarks/pseudoknot_benchmark_set.py at the top of the repository),
  evaluated as one strand, four strands, circular, and four strands with
  the third circular; benches/workloads/energy_workload.rs (an example
  target) for instruction counts
- RingDecomposition (for MultiPairTable): a structure as loops traced along
  the backbone (NestedLoop: its pairs in loop order and a Gap, closed, open
  from a 5' end to a 3' end, or through the wrap of a circle) and
  pseudoknots. All loops are traced, also those that meet their pairs out
  of order. Without crossing pairs, the loops are found pair by pair,
  walking between its two bases: a strand break ends an open loop, a strand
  wrap closes a circle. The loop graph (a vertex per loop, an edge per
  pair) decides what is strained: a Pseudoknot holds the
  loops of blocks (biconnected components) that share loops, with each
  segment strained if two visits of the same block enclose it, and the
  interior loops of its strained stems. EnergyModel::energy_of_structure,
  energy_of_nested_loop and energy_of_pseudoknot evaluate them
- Pseudoknots: evaluated with PKParams (defaults from a preliminary
  training on the HotKnots data sets); Fuzzyfold::with_pseudoknot_params; a warning
  when they are used at a temperature other than 37 C. Helix ends of the
  loops of a pseudoknot get their terminal penalty, and a mismatch
  (strained: mismatch_multi) or dangles
- Multi-stranded structures with pseudoknots: energies independent of the
  order of the strands. Fuzzyfold counts the duplex initiation per complex,
  one for each strand joined (strands minus complexes): a strand may be
  joined by a pseudoknot alone. The complexes are the components of the
  loop graph: loops minus pairs plus the cycles of the pseudoknots
  (Pseudoknot::cycles)
- NestedLoop is also the loop type of the move sets: split_loop,
  join_loop, extended, stretches (the unpaired bases between helix ends)
  and hairpin_size, also through the wrap of a circle
- Fuzzyfold::energy_terms, ViennaRNA::energy_terms: the energy of every
  loop and pseudoknot term (`ff-eval -v`); energy_of_structure evaluates
  only, it no longer logs the loops at log level Info
- check_sequence: does a structure fit a sequence (length, strand ends)?
  For input; evaluation assumes it without checking
- Circular strands in sequences: Base::SW (`^`), NucleotideVec::has_circular_strand
- Circular strands, evaluated alone as in ViennaRNA 2.6 (RNAeval -c):
  the loop through the wrap (Gap::Wrap) is the open circle (0), a hairpin,
  an interior or a multiloop. EnergyError::CircularStrand for a circular strand with other
  strands
- Circular strands with pseudoknots (Fuzzyfold): the loop through the wrap
  is closed by the backbone, strained once it lies on a cycle; a stem
  continues through the wrap (NestedLoop::is_stack); energies independent
  of where the circle is cut
- A circular strand together with linear strands (Fuzzyfold): the walk
  continues through the wrap at the first base of the circle; a circle
  joined to another strand lies on a cycle of the loop graph (its ring), so
  it is part of a pseudoknot. Several circles, also joined to each other:
  a segment may pass the wrap of another circle than the one of its loop
  (NestedLoop::segments copies the segments that pass a wrap). Energies
  independent of the order of the strands and of where each circle is cut
- Fuzzyfold: the open circle of n bases costs the initiation of a hairpin
  loop with n - 2 unpaired bases (the ring it closes), relative to the
  open chain; ViennaRNA keeps 0, as ViennaRNA 2.6

### Changed
- Examples of ViennaRNA (parameter sets, temperature) and Fuzzyfold
  (strands in either order).
- pastey instead of paste (unmaintained, RUSTSEC-2024-0436).
- The regression tests use the benchmark sets: tests/pseudoknot_energies.rs
  compares the energies of the pseudoknotted sets (lengths 50 and 100, all
  four variants, also the circles under rotation) with those recorded by
  benchmarks/fuzzyfold_energies.py, and tests/reference_energies.rs those of
  the nested sets (Turner 2004, Andronescu 2007, Mathews 2004) with
  ViennaRNA (benchmarks/viennarna_energies.py). tests/data is removed; the
  long hairpin, also at 25 C, is a unit test of ViennaRNA.
- The Display of loops and pseudoknots is plain text (no terminal colors);
  ff_energy no longer depends on colored, env_logger, itertools and
  rustc-hash
- NucleotideVec keeps its bases private (a sequence is checked: every
  strand has a base); TryFrom<Vec<Base>> builds one from bases, e.g. a
  modified copy. EnergyModel::energy_of_structure has no default
- ViennaRNA is the model of ViennaRNA 2.6 and does not evaluate pseudoknots
  (EnergyError::UnsupportedPseudoknot); ViennaRNA::default() uses Turner 2004.
  Fuzzyfold is the model under development (ViennaRNA loops, pseudoknots),
  Fuzzyfold::default() uses RNA_EXTENDED
- Indices are usize (NestedLoop), not NAIDX
- Base::I is Base::IN (modified bases have two-letter names); I is not a
  DNA base
- The module parameters::rna is parameters::rna_extended; K0 and
  parameters::T_REF are no longer public

### Removed
- LoopDecomposition and Decomposed (replaced by RingDecomposition), also
  for PairTable (PairTable is removed from ff_structure)
- NearestNeighborLoop, the loops of the old decomposition, with
  EnergyModel::energy_of_loop: NestedLoop is the loop type of energies and
  move sets. Joining an exterior and a joint exterior loop (a
  dissociation, NearestNeighborLoop::Disconnected) is left to future
  dissociation moves

## 0.5.0 - 2026-10-06
### Fixed
- GU stacking enthalpies.

### Added
- Inosine parameters.

### Changed
- Handling of closing penalties.

## [0.4.1] - 2026-06-18
### Fixed
- Revised duplicate hairpin free-energy parameters.

### Added
- is_stack for NearestNeighborLoop (helper for shift-classification).

### Changed
- Added Send + Sync supertrait for Energy-model.

## [0.4.0] - 2026-03-04
### Fixed
- Mayor compile-time boost.

### Changed
- Switched to supplying parameters in source / no more parameter files.

### Added
- Incorporated pseudouridine parameters (P).
- Added distinctions between RNA vs DNA, Thermodynamic vs Fitted parameters.
- Added fallback modes for smaller parameter-sets.

## [0.3.1] - 2026-02-06
### Fixed
- fix the evaluation of exterior-loops of length 1: '.' 

## [0.3.0] - 2026-01-13
### Added
- multifold evaluation (changes across many files)

### Changed
- Using FxHashMap instead of AHashMap for performance.
- Cargo bench updates.

### Removed
- Preliminary implementations for coaxial stacking.
