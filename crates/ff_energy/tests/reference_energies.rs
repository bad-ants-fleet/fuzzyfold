//! The ViennaRNA model against RNAeval: the parameter sets of Turner 2004
//! (also at another temperature), Andronescu 2007 and Mathews 2004 (DNA),
//! on structures of the benchmark set and a hairpin longer than the tables.
//! Away from 37 C only the hairpin: the model rounds rescaled parameters,
//! where ViennaRNA truncates them (energies differ by a few tenths of a
//! kcal/mol, see `ViennaRNA`).

use ff_structure::MultiPairTable;
use ff_energy::{EnergyModel, NucleotideVec, ViennaRNA};
use ff_energy::parameters::{DNA_MATHEWS_2004, RNA_ANDRONESCU_2007, RNA_TURNER_2004};

#[test]
fn energies_as_rnaeval() {
    let data = include_str!("data/reference_energies.txt");
    let mut tested = 0;
    for line in data.lines().filter(|l| !l.starts_with('#')) {
        let [params, celsius, seq, db, energy] = line.split_whitespace().collect::<Vec<_>>()[..] else {
            panic!("Malformed line: {line}");
        };
        let celsius: f64 = celsius.parse().unwrap();
        let (model, sequence) = match params {
            "turner2004" => (ViennaRNA::from_thermo_params(&RNA_TURNER_2004, celsius), NucleotideVec::try_from_rna(seq)),
            "andronescu2007" => (ViennaRNA::from_andronescu_params(&RNA_ANDRONESCU_2007), NucleotideVec::try_from_rna(seq)),
            "mathews2004" => (ViennaRNA::from_thermo_params(&DNA_MATHEWS_2004, celsius), NucleotideVec::try_from_dna(seq)),
            _ => panic!("Unknown parameters: {params}"),
        };
        let structure = MultiPairTable::try_from(db).unwrap();
        let e = model.energy_of_structure(&sequence.unwrap(), &structure).unwrap();
        assert_eq!(e, energy.parse::<i32>().unwrap(), "{params} {celsius} {seq} {db}");
        tested += 1;
    }
    assert_eq!(tested, 19);
}
