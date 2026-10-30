# Changelog ff_structure

All notable changes to this crate will be documented in this file.

## [development] - 0.6.0
### Added
- A README of the crate, with an example tested with the doctests.
- Pseudoknots: brackets `()[]{}<>` and `Aa`..`Zz`, 30 levels
- Circular strands: strand wrap `^` (DotBracket::STRAND_WRAP, MultiStruct::StrandWrap)
- MultiPairTable::is_nested: no two pairs cross; a strand wrap counts as
  a pair of the first and the last base of its circle, so a pair from a
  circle to another strand crosses
- ClosedRegions (From a MultiPairTable): closed regions and density, a
  structure analysis (the energy evaluation does not use it) The density
  of a circular strand depends on where the circle is cut (documented,
  with a test).
- crate documentation with examples

### Changed
- DotBracket is one byte: OPEN, CLOSE, UNPAIRED, STRAND_BREAK, open(level), close(level)
- MultiPairTable: two bytes per position, pseudoknots, access with get/set
- Indices are usize (MultiStruct::Paired, get, set, ClosedRegions); NAIDX
  is private to MultiPairTable, which stores positions in two bytes
- MultiPairTable -> DotBracketVec is TryFrom
- Empty strands are an error (a leading or trailing `+`, `++`)
- StructureError: Unmatched{Open,Close} carry the bracket; UnmatchedMulti* replaced by TooManyPseudoknotLevels, TooLong

### Removed
- PairTable (use MultiPairTable), StrandPairTable, LoopTable, PairSet, PairList, Pair, P1KEY
- dependencies nohash-hasher, serde

## [0.4.0] - 2026-10-06
### Added
- PairList

## [0.3.1] - 2026-01-13
### Added
- NAIDX indexing for PairTable

## [0.3.0] - 2026-01-13
### Added/Changed
- multi_pair_table.rs: new structures for MultiPairTable & StrandPairTable


