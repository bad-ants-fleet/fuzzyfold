# Changelog ff_energy

All notable changes to this crate will be documented in this file.

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

