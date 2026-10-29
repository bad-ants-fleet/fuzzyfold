use log::{debug, info};
use colored::*;
use clap::Args;
use clap::ValueEnum;
use anyhow::{bail, Result};

use ff_energy::Base;
use ff_energy::EnergyError;
use ff_energy::EnergyModel;
use ff_energy::Fuzzyfold;
use ff_energy::NucleotideVec;
use ff_energy::ViennaRNA;
use ff_energy::Pseudoknot;
use ff_energy::NestedLoop;
use ff_energy::RingDecomposition;
use ff_energy::parameters::RNAThermoParams;
use ff_energy::parameters::RNA_EXTENDED;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::RNA_ANDRONESCU_2007;
use ff_energy::parameters::DNA_MATHEWS_2004;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ModelChoice {
    /// The model of fuzzyfold (with pseudoknots).
    Fuzzyfold,
    /// The model of ViennaRNA 2.6 (without pseudoknots).
    Viennarna,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum RnaParams {
    Turner2004ext,
    Turner2004,
    Andronescu2007,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DnaParams {
    Mathews2004,
}

/// Free energy evaluation parameters.
#[derive(Debug, Args)]
pub struct EnergyModelArguments {
    /// Energy model
    #[arg(long, value_enum, default_value = "fuzzyfold")]
    pub model: ModelChoice,

    /// Temperature in Celsius
    #[arg(long, default_value = "37.0", allow_hyphen_values = true)]
    pub celsius: f64,

    /// Built-in RNA parameter set (default: turner2004ext, or turner2004
    /// with --model viennarna).
    #[arg(long, value_enum,
        num_args = 0..=1,
        value_name = "RNA_PRESET",
        conflicts_with_all = ["dna"])]
    pub rna: Option<Option<RnaParams>>,

    /// Built-in DNA parameter set (default: mathews2004).
    #[arg(long, value_enum,
        num_args = 0..=1,
        value_name = "DNA_PRESET",
        conflicts_with_all = ["rna"])]
    pub dna: Option<Option<DnaParams>>,
}

/// The energy model chosen on the command line.
pub enum Model {
    Fuzzyfold(Fuzzyfold),
    ViennaRNA(ViennaRNA),
}

impl Model {
    /// The energy terms of a structure (see `Fuzzyfold::energy_terms`).
    pub fn energy_terms<T: RingDecomposition>(&self, sequence: &[Base], structure: &T
    ) -> Result<Vec<(usize, String, i32)>, EnergyError> {
        match self {
            Model::Fuzzyfold(m) => m.energy_terms(sequence, structure),
            Model::ViennaRNA(m) => m.energy_terms(sequence, structure),
        }
    }
}

/// Log the energy terms of a structure, indented by level (`ff-eval -v`).
pub fn log_terms(terms: &[(usize, String, i32)]) {
    for (level, term, en) in terms {
        let en = format!("{:>6.2}", *en as f64 / 100.);
        match level {
            0 => info!("{:<41} {}", format!("{}:", colored_type(term)), en.green()),
            _ => info!("{:w$}{:<x$} {}", "", colored_type(term), en, w = 2 * level + 2, x = 39 - 2 * level),
        }
    }
}

/// A loop or pseudoknot description with its type (the first word) in color.
fn colored_type(term: &str) -> String {
    let Some((word, rest)) = term.split_once(' ') else { return term.to_string() };
    let word = match word {
        "Hairpin" | "Interior" => word.cyan(),
        "Multibr." | "Exterior" | "Circular" => word.cyan().bold(),
        "Jxterior" => word.red().bold(),
        "Pseudol." => word.magenta(),
        "Pseudok." => word.magenta().bold(),
        _ => return term.to_string(),
    };
    format!("{word} {rest}")
}

impl EnergyModelArguments {
    /// Warn that the ViennaRNA model is off for circular strands: it scores
    /// the open circle as 0, but charges the closure of every loop of a
    /// structured circle, relative to the open chain.
    pub fn warn_circular(&self, sequence: &NucleotideVec) {
        if matches!(self.model, ModelChoice::Viennarna) && sequence.has_circular_strand() {
            eprintln!("{} --model viennarna scores the open circle as 0 (as ViennaRNA 2.6), \
                which is not the reference state of the structured circles. Their \
                energies are relative to the open chain; use --model fuzzyfold.", 
                "Warning:".yellow().bold());
        }
    }

    pub fn build_model(&self) -> Result<Model> {
        debug!("{} {} C", "Celsius:".bold().red(), self.celsius);
        let default_rna = match self.model {
            ModelChoice::Fuzzyfold => RnaParams::Turner2004ext,
            ModelChoice::Viennarna => RnaParams::Turner2004,
        };
        let thermo: &'static RNAThermoParams = if let Some(dna_choice) = &self.dna {
            match dna_choice.unwrap_or(DnaParams::Mathews2004) {
                DnaParams::Mathews2004 => &DNA_MATHEWS_2004,
            }
        } else {
            match self.rna.flatten().unwrap_or(default_rna) {
                RnaParams::Turner2004 => &RNA_TURNER_2004,
                RnaParams::Turner2004ext => &RNA_EXTENDED,
                RnaParams::Andronescu2007 => {
                    if self.celsius != 37.0 {
                        bail!("Fitted parameters (andronescu2007) are given at 37 C only.");
                    }
                    return Ok(match self.model {
                        ModelChoice::Fuzzyfold => Model::Fuzzyfold(
                            Fuzzyfold::from_andronescu_params(&RNA_ANDRONESCU_2007)),
                        ModelChoice::Viennarna => Model::ViennaRNA(
                            ViennaRNA::from_andronescu_params(&RNA_ANDRONESCU_2007)),
                    });
                }
            }
        };
        Ok(match self.model {
            ModelChoice::Fuzzyfold => Model::Fuzzyfold(Fuzzyfold::from_thermo_params(thermo, self.celsius)),
            ModelChoice::Viennarna => Model::ViennaRNA(ViennaRNA::from_thermo_params(thermo, self.celsius)),
        })
    }
}

impl EnergyModel for Model {
    fn can_pair(&self, b1: Base, b2: Base) -> bool {
        match self {
            Model::Fuzzyfold(m) => m.can_pair(b1, b2),
            Model::ViennaRNA(m) => m.can_pair(b1, b2),
        }
    }

    fn min_hairpin_size(&self) -> usize {
        match self {
            Model::Fuzzyfold(m) => m.min_hairpin_size(),
            Model::ViennaRNA(m) => m.min_hairpin_size(),
        }
    }

    fn temperature(&self) -> f64 {
        match self {
            Model::Fuzzyfold(m) => m.temperature(),
            Model::ViennaRNA(m) => m.temperature(),
        }
    }

    fn energy_of_structure<T: RingDecomposition>(&self,
        sequence: &[Base],
        structure: &T
    ) -> Result<i32, EnergyError> {
        match self {
            Model::Fuzzyfold(m) => m.energy_of_structure(sequence, structure),
            Model::ViennaRNA(m) => m.energy_of_structure(sequence, structure),
        }
    }

    fn energy_of_pseudoknot(&self,
        sequence: &[Base],
        pseudoknot: &Pseudoknot
    ) -> Result<i32, EnergyError> {
        match self {
            Model::Fuzzyfold(m) => m.energy_of_pseudoknot(sequence, pseudoknot),
            Model::ViennaRNA(m) => m.energy_of_pseudoknot(sequence, pseudoknot),
        }
    }

    fn energy_of_nested_loop(&self,
        sequence: &[Base],
        nested: &NestedLoop
    ) -> Result<i32, EnergyError> {
        match self {
            Model::Fuzzyfold(m) => m.energy_of_nested_loop(sequence, nested),
            Model::ViennaRNA(m) => m.energy_of_nested_loop(sequence, nested),
        }
    }
}
