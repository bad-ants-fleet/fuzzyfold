# Changelog fuzzyfold

All notable changes to this crate will be documented in this file.

## [0.4.1] - 2026-02-26
## Changed
- Updates for energy evaluation (compiled parameters only, P support)
- Updates for rate model selection (Arrhenius only)

## [0.4.0] - 2026-02-06
## Added 
- ff-trajectory silent mode.
- ff-trajectory co-transcriptional mode.
- three-way and four-way shift move support for ff-explore and Metropolis rate model.

## Changed
- ff-trajectory / ff-timecourse input formats!
  Both programs now require an input structure. (A starting structure can now
  only be omitted in co-transcriptional mode.)
- ff-explore instead of ff-locmin
  Renamed to avoid confusion with different software. Supports shift moves.
- new default parameters in kinetic parsers.


## [0.3.0] - 2026-01-13
### Added
- ff-eval now supports multi-fold
- ff-locmin release
- Kawasaki support for ff-trajectory / ff-timecourse

### Changed
- ff-trajectory output formatting
- ff-timecourse k0 default for consistency

