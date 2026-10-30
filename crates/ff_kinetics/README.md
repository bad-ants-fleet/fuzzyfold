# ff_kinetics

Stochastic folding kinetics of nucleic acid secondary structures for
[fuzzyfold](https://github.com/bad-ants-fleet/fuzzyfold): the stochastic
simulation algorithm (Gillespie) on a move set of base pair additions and
deletions, three-way and four-way shift moves and A-to-I edits;
co-transcriptional folding; macrostates; time courses (`.tln`, `.nxy`)
and ensembles (`.drf`) of many trajectories; and neighborhoods of
structures.

```bash
cargo add ff_structure ff_energy ff_kinetics rand
```

```rust
use std::sync::Arc;
use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec};
use ff_kinetics::{Arrhenius, LoopNeighbors, MoveSet, SSA};

let model = Arc::new(Fuzzyfold::default());
let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
let start = MultiPairTable::try_from("...........").unwrap();

// Adding and deleting pairs at k0 = 1e5 per second; shift moves and edits
// are added with .with_three_way(k3ws), .with_four_way(k4ws) and
// .with_a_to_i(k_atoi).
let rates = Arrhenius::new(model.temperature(), 1e5);
let moves = LoopNeighbors::try_from((seq, &start, model, rates.move_types())).unwrap();
let mut ssa = SSA::from((moves, rates));

// Simulate 1 ms: time, structure and energy (kcal/mol) of each step.
ssa.simulate(&mut rand::rng(), 1e-3, |t, _waiting_time, _flux, m| {
    println!("{t:.2e} {} {:6.2}", m.current_structure(), m.current_energy() as f64 / 100.0);
    true
});
```

Kinetics support one strand without pseudoknots, linear or circular. See
the [documentation](https://docs.rs/ff_kinetics) for output times,
macrostates and time courses, and the
[workspace](https://github.com/bad-ants-fleet/fuzzyfold) for the command
line tools, the Python bindings and the citation.
