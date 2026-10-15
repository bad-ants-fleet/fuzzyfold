//! Standalone loop for profiling pseudo_findpath with samply/instruments.
//! Not a user-facing tool; lives in src/bin only during development.
//!
//! Usage: cargo build -p ff_findpath --bin ff-pseudo-findpath-bench-runner --profile bench
//!        samply record ./target/release/ff-pseudo-findpath-bench-runner

use ff_findpath::pseudo_findpath::findpath_pseudo;
use ff_energy::{parameters::{RNA_DP09, RNA_MT09}, ViennaRNA};
use std::hint::black_box;

const SEQ: &str = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCCCCCCC";
const TGT: &str = "(((((((((((((((((((([[[[[[[[[[[[[[[[[[[[))))))))))))))))))))..........]]]]]]]]]]]]]]]]]]]]";
const BEAM: usize = 50;
const REPS: usize = 50;

fn main() {
    let model = ViennaRNA::from_andrunescu_params(&RNA_MT09)
        .with_pseudoknot_params(RNA_DP09);

    for _ in 0..REPS {
        let result = findpath_pseudo(
            black_box(&model),
            black_box(SEQ),
            black_box(None),
            black_box(TGT),
            black_box(BEAM),
            black_box(None),
            false,
        ).unwrap();
        black_box(result);
    }
}
