use log::debug;
use colored::*;
use clap::Args;
use clap::ValueEnum;

use ff_energy::ViennaRNA;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::RNA_TURNER_2004_EXTENDED;
use ff_energy::parameters::RNA_ANDRONESCU_2007;
use ff_energy::parameters::DNA_MATHEWS_2004;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum RnaParams {
    Turner2004,
    Turner2004ext,
    Andronescu2007,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum DnaParams {
    Mathews2004,
}

/// Free energy evaluation parameters.
#[derive(Debug, Args)]
pub struct EnergyModelArguments {
    /// Temperature in Celsius
    #[arg(short, long, default_value = "37.0")]
    pub temperature: f64,

    /// Built-in RNA parameter set (default: turner2004).
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

impl EnergyModelArguments {
    pub fn build_model(&self) -> ViennaRNA {
        debug!("{} {} °C", "Temperature:".bold().red(), self.temperature);
        if let Some(rna_choice) = &self.rna {
            let preset = rna_choice.unwrap_or(RnaParams::Turner2004);
            match preset {
                RnaParams::Turner2004 => {
                    ViennaRNA::from_thermo_params(&RNA_TURNER_2004, self.temperature)
                },
                RnaParams::Turner2004ext => {
                    ViennaRNA::from_thermo_params(&RNA_TURNER_2004_EXTENDED, self.temperature)
                },
                RnaParams::Andronescu2007 => {
                    if self.temperature != 37.0 {
                        panic!("Cannot change temperature for fitted parameters!");
                    }
                    ViennaRNA::from_fitted_params(&RNA_ANDRONESCU_2007)
                },
            }
        } else if let Some(dna_choice) = &self.dna {
            let preset = dna_choice.unwrap_or(DnaParams::Mathews2004);
            match preset {
                DnaParams::Mathews2004 => {
                    ViennaRNA::from_thermo_params(&DNA_MATHEWS_2004, self.temperature)
                },
            }
        } else {
            ViennaRNA::from_thermo_params(&RNA_TURNER_2004, self.temperature)
        }
    }
}


