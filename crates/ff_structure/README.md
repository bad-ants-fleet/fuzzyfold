# ff_structure

Nucleic acid secondary structures for [fuzzyfold](https://github.com/bad-ants-fleet/fuzzyfold):
dot-bracket strings with pseudoknots and several strands (`DotBracketVec`),
pair tables (`MultiPairTable`), and the closed regions and pseudoknot
density of a structure (`ClosedRegions`).

```bash
cargo add ff_structure
```

```rust
use ff_structure::{ClosedRegions, DotBracketVec, MultiPairTable, MultiStruct};

// Two strands (+), with a pseudoknot ([]).
let db = DotBracketVec::try_from("((.[[+))..]]").unwrap();
let pairs = MultiPairTable::try_from(&db).unwrap();
assert_eq!(pairs.get(0), MultiStruct::Paired(7));
assert_eq!(pairs.get(5), MultiStruct::StrandBreak);
assert!(!pairs.is_nested());
assert_eq!(ClosedRegions::from(&pairs).density(), 2);
```

Indices are positions in the structure string: a strand break (`+`) or
strand wrap (`^`, after a circular strand) takes one position. See the
[documentation](https://docs.rs/ff_structure) and the
[workspace](https://github.com/bad-ants-fleet/fuzzyfold) for the other
crates, the command line tools and the citation.
