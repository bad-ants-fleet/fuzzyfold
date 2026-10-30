use clap::Args;
use anyhow::Result;
use ff_kinetics::Arrhenius;
use ff_kinetics::OutputTimes;
use ff_kinetics::extension_schedule;

/// Rate model parameter parsing.
#[derive(Debug, Args)]
pub struct RateModelArguments {
    /// Rate constant for add/delete moves.
    #[arg(long, default_value_t = 1e5)]
    pub k0: f64,

    /// Rate constant for three-way shift moves (optional, default = off).
    #[arg(long)]
    pub k3ws: Option<f64>,

    /// Rate constant for four-way shift moves (optional, default = off).
    #[arg(long)]
    pub k4ws: Option<f64>,
}

impl RateModelArguments {
    /// Validate that all parameters make sense.
    pub fn build_model(&self, celsius: f64) -> Arrhenius {
        Arrhenius::new(celsius, self.k0, self.k3ws, self.k4ws)
    }
}

/// The simulation times: extensions (co-transcriptional) and output times.
///
/// The output times are linear up to --t-sep (--t-lin times) and
/// logarithmic after it (--t-log times), up to the end of the simulation.
///  - full length: the simulation ends at --t-end; --t-sep defaults to
///    1e-4 s.
///  - co-transcriptional (--t-ext): one extension every --t-ext, and the
///    simulation ends --t-end after the full length. --t-sep defaults to
///    the end of transcription, with one output time per extension.
#[derive(Debug, Args)]
pub struct TimeArguments {
    /// Time per extension during transcription (co-transcriptional
    /// simulation).
    #[arg(long)]
    pub t_ext: Option<f64>,

    /// Simulation time at the full-length sequence.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// End of the linear output times. Default: 1e-4, or the end of
    /// transcription.
    #[arg(long)]
    pub t_sep: Option<f64>,

    /// Number of linear output times. Default: 50, or one per extension.
    #[arg(long)]
    pub t_lin: Option<usize>,

    /// Number of logarithmic output times.
    #[arg(long, default_value_t = 50)]
    pub t_log: usize,
}

impl TimeArguments {
    /// The output times of a simulation with num_ext extensions.
    pub fn output_times(&self, num_ext: usize) -> Result<OutputTimes> {
        if self.t_ext.is_some() && num_ext == 0 {
            anyhow::bail!("--t-ext is given, but the start structure has the length of the sequence");
        }
        let transcription = self.t_ext.map_or(0.0, |t| t * num_ext as f64);
        let (t_sep, t_lin) = match (self.t_ext, self.t_sep) {
            (Some(_), None) => (transcription, self.t_lin.unwrap_or(num_ext)),
            (_, t_sep) => (t_sep.unwrap_or(1e-4), self.t_lin.unwrap_or(50)),
        };
        OutputTimes::lin_log(t_sep, transcription + self.t_end, t_lin, self.t_log)
            .map_err(anyhow::Error::msg)
    }

    /// The times of `SSA::co_simulate` up to the end of the output times.
    pub fn schedule(&self, num_ext: usize, output: &OutputTimes) -> Result<Vec<f64>> {
        extension_schedule(self.t_ext, num_ext, output.end()).map_err(anyhow::Error::msg)
    }
}
