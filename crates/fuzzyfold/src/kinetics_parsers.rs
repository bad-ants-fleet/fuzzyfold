use clap::Args;
use anyhow::bail;
use anyhow::Result;
use ff_kinetics::Arrhenius;
use plotters::prelude::*;

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

/// Timeline parameter parsing.
///
/// Any timeline is composed of two regimes: 
///  - first a linear regime [0 .. t-seq] 
///  - second a logarithmic regime [t-sep .. t-end].
///
/// Each regime is divided into equally spaced time-points for analysis, which
/// are specified by parameters --t-lin (>=1) and --t-log (>=0) respectively.
///
/// We distinguish two modes:
///  - full-length: When the sequence length remains constant, then --t-sep
///    defaults to 1/k0, i.e. the expected waiting time for the fastest reactions.
///    All other parameters have constant defaults and t-end is set explicitly
///    with --t-end.
///
///  - co-transcriptional: When sequence length changes over time, then the
///    additional parameter --t-ext sets the simulation time at each nucleotide.
///    In this mode, --t-sep defaults to the end-of-transcription, and 
///    t-end = #extensions * --t-ext + --t-end -- that is, the --t-end parameter
///    only effects the simulation time at the last nucleotide.
///
#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// Extension time during transcription. Switches the default full-length to
    /// co-transcriptional timeline mode.
    #[arg(long)]
    pub t_ext: Option<f64>,   

    /// Simulation time at the full-length sequence.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Sets the end of a linear output time regime and start of a logarithmic one, s.t. 0 < t-sep
    /// <= t-end. In full-length mode, defaults to 1/k0. In a co-transcriptional simulation,
    /// defaults to the end of transcription.
    #[arg(long)]
    pub t_sep: Option<f64>,   

    /// Number of time points on the linear scale. Default = 0 is a hack to default to 1 in full-length mode,
    /// and end-of-transcription in co-transcriptional mode.
    #[arg(long)]
    pub t_lin: Option<usize>, 

    /// Number of time points on the logarithmic scale.
    #[arg(long, default_value_t = 50)]
    pub t_log: usize,
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    /// Either: t-ext is none: classic-mode
    pub fn validate(&mut self, k0: f64, num_ext: usize) -> Result<()> {
        if (num_ext == 0) != self.t_ext.is_none() {
            // needs better bail warning for different cases.
            bail!("Inconsistent input, if t_ext is set, then there must be extensions. If there are extensions, t_ext must be set!");
        }

        // Set default values for t_sep in case it is not set by user.
        if self.t_lin.is_none() {
            // full-length mode or if t-sep is set.
            if self.t_ext.is_none() || self.t_sep.is_some() { 
                self.t_lin = Some(50);
            } else { // co-transcriptional mode
                self.t_lin = Some(num_ext);
            }
        } 

        let total_ext = self.t_ext.map(|t| t * num_ext.as_f64());
        match (self.t_sep, total_ext) {
            // Defaults if t_sep is not set by the user.
            (None, None) => self.t_sep = Some(10.0 / k0),  // full-length mode
            (None, Some(ext)) => self.t_sep = Some(ext),   // co-transcriptional mode
            // Verify user-set t_sep.
            (Some(sep), None) => {
                if self.t_end <= sep {
                    bail!("t_end ({}) must be greater than t_sep ({})", self.t_end, sep);
                }
            }
            (Some(sep), Some(ext)) => {
                if ext + self.t_end <= sep {
                    bail!("t_sep ({}) must be smaller than the total simulation time ({})", sep, ext + self.t_end);
                }
            }
        }
 
        Ok(())
    }

    pub fn get_output_times(&self, num_ext: usize) -> Result<Vec<f64>> {
        let t_end = self.t_end;
        let t_lin = self.t_lin.expect("t-lin has to be set during validation!");
        let t_log = self.t_log;
        let t_sep = self.t_sep.expect("t-sep has to be set during validation!");
        let end = if let Some(t_ext) = self.t_ext {
            t_ext * num_ext.as_f64() + t_end
        } else { t_end };

        if t_lin == 0 {
            if t_log != 1 {
                bail!("If t_lin == 0, then t_log = 1! (A special hack to support single-output mode).");
            }
            return Ok(vec![end]);
        }

        let mut times = vec![0.0];
        let start = *times.last().unwrap();
        let step = t_sep / t_lin as f64;
        for i in 1..=t_lin {
            times.push(start + i as f64 * step);
        }

        // Logarithmic tail: append 't_log logarithmic timepoints between t-sep...t_end
        let start = *times.last().unwrap();
        let log_start = start.ln();
        let log_end = end.ln();
        for i in 1..t_log {
            let frac = i as f64 / t_log as f64;
            let value = (log_start + frac * (log_end - log_start)).exp();
            times.push(value);
        }
        times.push(end);

        Ok(times)
    }
}


