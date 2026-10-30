//! The ViennaRNA model against ViennaRNA itself: the energies of the nested
//! benchmark sets (lengths 50 and 100) with the parameter sets Turner 2004,
//! Andronescu 2007 and Mathews 2004 (DNA), at 37 C, as recorded by
//! `benchmarks/viennarna_energies.py` (`benches/data/*_viennarna.txt`).

use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, NucleotideVec, ViennaRNA};
use ff_energy::parameters::{DNA_MATHEWS_2004, RNA_ANDRONESCU_2007, RNA_TURNER_2004};

const SETS: [(&str, &str); 2] = [
    (include_str!("../benches/data/benchmark_random_structures_len50.vrna"),
     include_str!("../benches/data/benchmark_random_structures_len50_viennarna.txt")),
    (include_str!("../benches/data/benchmark_random_structures_len100.vrna"),
     include_str!("../benches/data/benchmark_random_structures_len100_viennarna.txt")),
];

#[test]
fn energies_as_viennarna() {
    let turner = ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0);
    let andronescu = ViennaRNA::from_andronescu_params(&RNA_ANDRONESCU_2007);
    let mathews = ViennaRNA::from_thermo_params(&DNA_MATHEWS_2004, 37.0);
    let mut tested = 0;
    for (set, energies) in SETS {
        let records: Vec<&str> = set.lines().collect();
        let energies: Vec<&str> = energies.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(records.len(), 3 * energies.len());
        for (record, line) in records.chunks(3).zip(energies) {
            let (seq, db) = (record[1], record[2]);
            let expected: Vec<i32> = line.split_whitespace().map(|e| e.parse().unwrap()).collect();
            let mpt = MultiPairTable::try_from(db).unwrap();
            let rna = NucleotideVec::try_from_rna(seq).unwrap();
            let dna = NucleotideVec::try_from_dna(&seq.replace('U', "T")).unwrap();
            let found = [
                turner.energy_of_structure(&rna, &mpt).unwrap(),
                andronescu.energy_of_structure(&rna, &mpt).unwrap(),
                mathews.energy_of_structure(&dna, &mpt).unwrap(),
            ];
            assert_eq!(found[..], expected[..], "turner2004 andronescu2007 mathews2004\n{seq}\n{db}");
            tested += 1;
        }
    }
    assert_eq!(tested, 200);
}
