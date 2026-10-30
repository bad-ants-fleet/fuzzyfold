//! Stochastic folding kinetics of nucleic acid secondary structures:
//! simulations with the stochastic simulation algorithm (Gillespie),
//! co-transcriptional folding, macrostates, and neighborhoods of
//! structures.
//!
//! Three parts make a simulation:
//! - A [`MoveSet`] proposes the moves of the current structure with their
//!   energy changes, and applies them. [`LoopNeighbors`] is the move set of
//!   fuzzyfold: adding and deleting a pair, and optionally three-way and
//!   four-way shift moves ([`MoveTypes`]). It keeps the loops of the
//!   structure (its [`LoopTable`]) and updates only the loops a move
//!   changes. A move set that also extends the strand at its 3' end
//!   implements [`Elongation`]. With [`MoveTypes::a_to_i`], the move set
//!   also edits an A opposite a C in a 1x1 interior loop to I, which pairs
//!   with the C ([`Move::Edit`]): the sequence changes during a simulation
//!   ([`MoveSet::current_sequence`]).
//! - A [`RateModel`] gives the rate of a move. [`Arrhenius`] uses one rate
//!   constant per move type: the constant for a move that does not raise
//!   the energy, else the constant times exp(-dE / RT); for shift moves, dE
//!   is the barrier of the move. Shift moves and edits are off unless their
//!   rate constant is set ([`Arrhenius::with_three_way`] and others).
//! - [`SSA`] simulates: it keeps the rates of all moves in a tree, draws
//!   the waiting time and the next move, and applies the move.
//!
//! ```
//! use std::sync::Arc;
//! use ff_structure::MultiPairTable;
//! use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
//! use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, SSA};
//!
//! let model = Arc::new(Fuzzyfold::default());
//! let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
//! let start = MultiPairTable::try_from("...........").unwrap();
//!
//! // Adding and deleting pairs at k0 = 1e5 per second.
//! let rates = Arrhenius::new(model.temperature(), 1e5);
//! let moves = LoopNeighbors::try_from((seq, &start, model, rates.move_types())).unwrap();
//! let mut ssa = SSA::from((moves, rates));
//!
//! // Simulate 1 ms. The callback gets each structure with the time it is
//! // reached; returning false would stop the simulation.
//! ssa.simulate(&mut rand::rng(), 1e-3, |t, _waiting_time, _flux, m| {
//!     println!("{t:.2e} {} {:6.2}", m.current_structure(), m.current_energy() as f64 / 100.0);
//!     true
//! });
//! ```
//!
//! # Units
//!
//! Energies are integers in units of 10 cal/mol, as in `ff_energy`; rates
//! are per second and times in seconds.
//!
//! # Co-transcriptional folding
//!
//! A simulation may start from a structure shorter than the sequence: the
//! strand grows by one unpaired base at its 3' end every t_ext
//! ([`SSA::co_simulate`], with the times from [`extension_schedule`]).
//! Circular strands are simulated at their full length.
//!
//! # Output times, macrostates and time courses
//!
//! [`SSA::sample`] records the state at given [`OutputTimes`]: any
//! increasing list, linear and then logarithmic ([`OutputTimes::lin_log`]),
//! one per extension ([`OutputTimes::transcription`]), or the default.
//! A [`MacrostateRegistry`] holds sets of structures (macrostates, read
//! from files) and classifies a structure; a [`timecourse::Timecourse`]
//! counts the macrostates of many trajectories at the output times and
//! writes them as `.tln` and `.nxy` files; an [`ensemble::Ensemble`]
//! counts their structures and writes them as a `.drf` file.
//!
//! ```
//! use std::sync::Arc;
//! use ff_structure::MultiPairTable;
//! use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
//! use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, SSA};
//! use ff_kinetics::{MacrostateRegistry, OutputTimes, extension_schedule};
//! use ff_kinetics::timecourse::Timecourse;
//!
//! let model = Arc::new(Fuzzyfold::default());
//! let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
//! let start = MultiPairTable::try_from("...........").unwrap();
//! let rates = Arrhenius::new(model.temperature(), 1e5);
//! let moves = LoopNeighbors::try_from((seq.clone(), &start, model.clone(), rates.move_types())).unwrap();
//!
//! // A macrostate in the format of a file: name, sequence, structures.
//! // Structures in no macrostate are "Unassigned" (index 0).
//! let mut registry = MacrostateRegistry::from((seq, model));
//! let hairpin = ">hairpin\nGGGGAAACCCC\n((((...))))\n.(((...))).\n";
//! registry.insert_from_reader(hairpin.as_bytes(), "hairpin", 0).unwrap();
//!
//! // 50 times evenly up to 0.1 ms, then 50 evenly in log space up to 1 s.
//! let output = OutputTimes::default();
//! // At full length from the start: no extensions, one run to the end.
//! let schedule = extension_schedule(None, 0, output.end()).unwrap();
//!
//! let mut timecourse = Timecourse::new(output.times(), registry.len());
//! for _ in 0..10 {
//!     let mut ssa = SSA::from((moves.clone(), rates));
//!     let trajectory = ssa.sample(&mut rand::rng(), &schedule, &output,
//!         |m| registry.classify(m.current_sequence(), &m.current_structure()));
//!     timecourse.add_trajectory(&trajectory);
//! }
//!
//! // The number of trajectories in each macrostate at the end.
//! let (names, counts) = timecourse.counts_by_name(&registry);
//! for (name, count) in names.iter().zip(counts.last().unwrap()) {
//!     println!("{name} {count}");
//! }
//! ```
//!
//! # Neighborhoods
//!
//! [`enum_neighbors::enumerate_neighbors`] lists the structures reachable
//! from the current structure of a move set within an energy range, and
//! optionally a number of moves, for example to define macrostates around
//! local minima.
//!
//! # Restrictions and errors
//!
//! Kinetics support a single strand without pseudoknots, linear or
//! circular (at full length); building a move set for anything else fails
//! with [`MoveSetError`].

pub mod timecourse;
pub mod ensemble;
pub mod enum_neighbors;

mod rate_tree;
mod rate_model;
mod stochastic_simulation;
mod macrostates;
mod movesets;

pub use rate_model::*;
pub use stochastic_simulation::*;
pub use macrostates::*;
pub use movesets::*;

/// The example of the README, tested with the doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExample;
