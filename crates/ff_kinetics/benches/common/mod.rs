//! Inputs shared by the benchmarks.

use std::io::BufRead;
use ff_structure::MultiPairTable;
use ff_energy::NucleotideVec;

/// The path of a file in benches/data.
pub fn data(file: &str) -> String {
    format!("{}/benches/data/{}", env!("CARGO_MANIFEST_DIR"), file)
}

/// The sequences and structures of benchmark_random_structures_len{len}.vrna.
pub fn load_inputs(len: usize) -> Vec<(NucleotideVec, MultiPairTable)> {
    let path = data(&format!("benchmark_random_structures_len{len}.vrna"));
    let file = std::fs::File::open(&path).expect("Cannot open input file");
    let lines: Vec<String> = std::io::BufReader::new(file).lines()
        .map(|l| l.unwrap())
        .collect();
    lines.chunks(3).map(|c| {
        assert!(c[0].starts_with('>'), "Malformed benchmarking input.");
        (NucleotideVec::try_from(c[1].as_str()).unwrap(),
         MultiPairTable::try_from(c[2].as_str()).unwrap())
    }).collect()
}
