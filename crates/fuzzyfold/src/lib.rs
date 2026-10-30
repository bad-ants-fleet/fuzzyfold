//! # fuzzyfold
//!
//! The main entry point for the fuzzyfold nucleic acid folding package, with
//! interfaces to multiple crates that are organized as a workspace. At this
//! level, you can also find argument parsers for the binaries shipped with the
//! fuzzyfold package.
//!
//! This crate re-exports the library crates as `fuzzyfold::structure`,
//! `fuzzyfold::energy` and `fuzzyfold::kinetics`; the binaries use them
//! through these paths. Programs that only need the libraries can depend
//! on ff_structure, ff_energy and ff_kinetics directly.
//!
//! # Command line tools
//!
//! Installed with `cargo install fuzzyfold`. Each tool reads a sequence
//! (and a structure) from a file or stdin, and `--help` lists its options.
//! Energies are in kcal/mol, times in seconds with the default rate
//! constant `--k0 1e5`. The commands below run in
//! [`examples/ff-timecourse`](https://github.com/bad-ants-fleet/fuzzyfold/tree/development/examples/ff-timecourse);
//! the READMEs in `examples/` explain them and their output.
//!
//! | Tool | Does |
//! |---|---|
//! | `ff-eval` | the free energy of a structure, also with pseudoknots and several strands |
//! | `ff-randseq` | random sequences |
//! | `ff-trajectory` | one stochastic trajectory, also co-transcriptional (`--t-ext`) |
//! | `ff-ensemble` | the structures of many trajectories at the output times (`.drf`) |
//! | `ff-explore` | the structures within an energy range of a structure: a macrostate |
//! | `ff-timecourse` | the macrostate occupancies of many trajectories (`.tln`, `.nxy`, `.svg`) |
//!
//! ```text
//! # The energy of a pseudoknot; -v lists its loops and terms.
//! printf ">pk\nGCGAUUUCUGACCGCUUUUUUGUCAG\n[[[....(((((]]]......)))))\n" | ff-eval
//!
//! # A trajectory of a random 50-nt sequence, from the open chain, for 1 s.
//! ff-randseq -l 50 --eval | ff-trajectory --t-end 1
//!
//! # The structures of 100 trajectories from a start structure, as dld1_lm3.drf.
//! cat dld1_lm3.na | ff-ensemble -n 100 --t-end 1 --output dld1_lm3
//!
//! # A macrostate: the structures within 2.8 kcal/mol, sorted by energy.
//! cat dld1_lm3.na | ff-explore --delta 2.8 --sorted --rna turner2004 > lm3.ms
//!
//! # The occupancy of the macrostates dld1*.ms in 100 trajectories up to 1 s.
//! cat dld1_lm3.na | ff-timecourse --macrostates dld1*.ms --t-end 1 -n 100 --output dld1_lm3_t1
//! ```


/// ff_structure: representations of secondary structures.
pub use ff_structure as structure;

/// ff_energy: sequences, loop decompositions and energy models.
pub use ff_energy as energy;

/// ff_kinetics: move sets, rate models, stochastic simulation, macrostates
/// and time courses.
pub use ff_kinetics as kinetics;

/// The citation, shown at the end of `--help` of the command line tools.
pub const CITATION: &str = "If you use fuzzyfold, please cite: Stefan Badelt (2026): fuzzyfold: \
    a high-performance framework for stochastic RNA folding kinetics. bioRxiv, \
    https://doi.org/10.64898/2026.06.17.732885";

/// Let a closed pipe end a command line tool quietly, as it ends other
/// tools (`ff-trajectory | head`). Rust ignores SIGPIPE, so printing to a
/// closed pipe would panic. Call it first in main.
pub fn restore_sigpipe() {
    // SAFETY: sets a signal disposition, before any other thread exists.
    #[cfg(unix)]
    unsafe {
        libc::signal(libc::SIGPIPE, libc::SIG_DFL);
    }
}

/// Various flavors of handling sequence/structure input.
pub mod input_parsers;

/// Exposing the currently supported parameters of fuzzyfold's energy models.
pub mod energy_parsers;

/// Exposing the currently supported parameters of fuzzyfold's rate models and simulation parameters.
pub mod kinetics_parsers;

