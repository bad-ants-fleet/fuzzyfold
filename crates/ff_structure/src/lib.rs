mod error;
mod dotbracket;
mod pair_table;
mod multi_pair_table;
mod loop_table;
mod pair_set;

pub use error::*;
pub use dotbracket::*;
pub use pair_table::*;
pub use multi_pair_table::*;
pub use loop_table::*;
pub use pair_set::*;


/// Nucleic Acid INdeX: we use `u16` (0 to 65k), which is plenty for nucleic acids.
/// Should you ever want to fold longer sequences, beware that `P1KEY` needs to
/// be *twice as large* (in bits) as `NAIDX`, since pairs `(NAIDX, NAIDX)` are
/// compacted into one `P1KEY`.
pub type NAIDX = u16;

/// Pair key. Must be >= 2×`NAIDX` in bit width so we can safely pack two indices.
pub type P1KEY = u32;

/// Compile-time sanity check: 2×NAIDX bits must fit into P1KEY.
const _: () = {
    debug_assert!(2 * NAIDX::BITS <= P1KEY::BITS);
};


