# Changelog fuzzyfold

All notable changes to this crate will be documented in this file.

## [development]
## Added
- --k-AtoI: the rate of A-to-I editing in ff-trajectory, ff-timecourse and
  ff-ensemble (an A opposite a C in a 1x1 interior loop becomes I and pairs).
  ff-trajectory prints the sequence (as far as transcribed) when an edit
  changes it, with the energy and times of the new state; ff-ensemble
  writes the sequence as a column of the .drf file (a structure keeps its id
  as the strand grows, not when it is edited); a macrostate file of
  ff-timecourse may give an edited sequence (I for A)
- --model fuzzyfold|viennarna: the energy model (default fuzzyfold);
  with viennarna, the default RNA parameters are turner2004

- A warning for circular strands with --model viennarna: it scores the
  open circle as 0, which is not the reference state of structured circles
- Circular strands in ff-timecourse, ff-trajectory, ff-ensemble and
  ff-explore (at full length)

## Changed
- TimelineParameters is TimeArguments; it gives OutputTimes (ff_kinetics).
  --t-sep defaults to 1e-4 s at full length, no longer 10/k0 (the same for
  the default k0); --t-lin 0 (a single output time) is no longer
  supported
- ff-timecourse counts the trajectories as they finish, instead of keeping
  one timeline per simulation in memory
- andronescu2007 at a temperature other than 37 C is an error, not a panic
- Pseudoknot energies from the ring decomposition (see ff_energy): a loop
  enclosing a kiss is strained, no longer a regular multiloop; `ff-eval -v`
  lists the loops and, for a pseudoknot, its terms (energy_parsers::log_terms).
  Other tools no longer log energy terms with -v; the loop types are
  colored here, not in ff_energy

## [0.5.0] - 2026-10-06
## Added
- ff-ensemble to produce detailed output from multiple simulations in DRF format.
- co-transcriptional folding support for ff-timecourse.
- t-sep commandline parameter for timeline construction.
- t-lin/t-log logic adapted for timeline construction.

## Changed
- t-ext now explicit for co-transcriptional folding. 

## [0.4.2] - 2026-02-26
## Added
- fasta flag for ff-trajectory
- title option for ff-timecourse

## Changed
- ff-timecourse output format is now *.nxy
- Now allowing negative temperature (Celsius).

## [0.4.1] - 2026-02-26
## Changed
- Updates for energy evaluation (compiled parameters only, P support).
- Parameters default to "extended" with additional haipin energies and modifications.
- Updates for rate model selection (Arrhenius only).
- Restricted alphabets for RNA {A, C, G, U, P} and DNA {A, C, G, T}.

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

