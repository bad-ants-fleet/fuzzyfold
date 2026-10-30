# ff_energy

Nearest neighbor free energies of nucleic acid secondary structures for
[fuzzyfold](https://github.com/bad-ants-fleet/fuzzyfold): `ViennaRNA` (the
model of ViennaRNA 2.6, for comparisons) and `Fuzzyfold` (the model under
development: also pseudoknots, several strands independent of their order,
and circular strands), with the parameter sets Turner 2004, Andronescu 2007,
Mathews 2004 (DNA) and an extended RNA set.

```bash
cargo add ff_energy
```

```rust
use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, Fuzzyfold, NucleotideVec, ViennaRNA};

let seq = NucleotideVec::try_from_rna("GGGGAAACCCC").unwrap();
let hairpin = MultiPairTable::try_from("((((...))))").unwrap();

// Energies are in units of 10 cal/mol: -450 is -4.50 kcal/mol.
assert_eq!(ViennaRNA::default().energy_of_structure(&seq, &hairpin).unwrap(), -450);

// A pseudoknot, with the extended RNA parameters at 37 C.
let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
let pk = MultiPairTable::try_from("[[[....(((((]]]......)))))").unwrap();
println!("{}", Fuzzyfold::default().energy_of_structure(&seq, &pk).unwrap());
```

See the [documentation](https://docs.rs/ff_energy) and the
[workspace](https://github.com/bad-ants-fleet/fuzzyfold) for the other
crates, the command line tools and the citation.
