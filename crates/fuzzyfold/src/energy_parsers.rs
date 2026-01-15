use log::debug;
use colored::*;
use clap::Args;
use clap::ValueEnum;
use std::path::PathBuf;

use ff_energy::ViennaRNA;
use ff_energy::parameters::RNA_TURNER_2004;
use ff_energy::parameters::RNA_TURNER_2004_EXT;
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
        conflicts_with_all = ["dna", "na_params"])]
    pub rna: Option<Option<RnaParams>>,

    /// Built-in DNA parameter set (default: mathews2004).
    #[arg(long, value_enum, 
        num_args = 0..=1, 
        value_name = "DNA_PRESET",
        conflicts_with_all = ["rna", "na_params"])]
    pub dna: Option<Option<DnaParams>>,

    /// Explicit parameter file (overrides RNA/DNA presets).
    #[arg(long, value_name = "FILE",
        conflicts_with_all = ["rna", "dna"])]
    pub na_params: Option<PathBuf>,
}

impl EnergyModelArguments {
    pub fn build_model(&self) -> ViennaRNA {
        debug!("{} {} °C", "Temperature:".bold().red(), self.temperature);

        let mut model = if let Some(path) = &self.na_params {
            debug!("{} {:?}", "Polymer:".bold().red(), path);
            ViennaRNA::from_parameter_file(path)
                .expect("Failed to load parameter file")
        } else if let Some(rna_choice) = &self.rna {
            // RNA mode
            let preset = rna_choice.unwrap_or(RnaParams::Turner2004);
            let (data, name) = match preset {
                RnaParams::Turner2004 => (RNA_TURNER_2004, "rna_turner_2004"),
                RnaParams::Turner2004ext => (RNA_TURNER_2004_EXT, "rna_turner_2004_ext"),
                RnaParams::Andronescu2007 => (RNA_ANDRONESCU_2007, "rna_andronescu_2007"),
            };
            debug!("{} {}", "Polymer:".bold().red(), name);
            ViennaRNA::from_parameter_str(data).unwrap()

        } else if let Some(dna_choice) = &self.dna {
            // DNA mode
            let preset = dna_choice.unwrap_or(DnaParams::Mathews2004);
            let (data, name) = match preset {
                DnaParams::Mathews2004 => (DNA_MATHEWS_2004, "dna_mathews_2004"),
            };
            debug!("{} {}", "Polymer:".bold().red(), name);
            ViennaRNA::from_parameter_str(data).unwrap()
        } else {
            debug!("{} rna_turner_2004", "Polymer:".bold().red());
            ViennaRNA::default()
        };

        model.set_temperature(self.temperature);
        model
    }
}


