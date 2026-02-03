use clap::Args;
use clap::ValueEnum;
use anyhow::bail;
use anyhow::Result;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum RateModelKind {
    Metropolis,
    Kawasaki,
}

#[derive(Debug, Args)]
pub struct RateModelArguments {
    /// Rate constant for add/delete moves.
    #[arg(long, default_value_t = 1e6)]
    pub k0: f64,

    /// Rate constant for three-way shift moves (optional, default = off).
    #[arg(long)]
    pub k3ws: Option<f64>,

    /// Rate constant for four-way shift moves (optional, default = off).
    #[arg(long)]
    pub k4ws: Option<f64>,

    /// Select a rate model.
    #[arg(long, value_enum, default_value = "metropolis")]
    pub rate_model: RateModelKind,
}

#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// The last time point of the linear scale.
    #[arg(long, default_value_t = 1e-7)]
    pub t_ext: f64,

    /// Simulation stop time.
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Number of time points on the linear scale: [0..t-ext]
    #[arg(long, default_value_t = 10)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale: [t-ext..t-end]
    #[arg(long, default_value_t = 100)]
    pub t_log: usize,
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        if self.t_end <= self.t_ext {
            bail!("t_end ({}) must be greater than t_ext ({})", self.t_end, self.t_ext);
        }
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        Ok(())
    }

    pub fn get_output_times(&self) -> Vec<f64> {
        let t_end = self.t_end;
        let t_ext = self.t_ext;
        let t_lin = self.t_lin;
        let t_log = self.t_log;
        let mut times = vec![0.0];

        // Linear segments: append `t_lin` evenly spaced points
        let start = *times.last().unwrap();
        let step = t_ext / t_lin as f64;
        for i in 1..=t_lin {
            times.push(start + i as f64 * step);
        }

        // Logarithmic tail
        let start = *times.last().unwrap();
        let log_start = start.ln();
        let log_end = t_end.ln();
        for i in 1..t_log {
            let frac = i as f64 / t_log as f64;
            let value = (log_start + frac * (log_end - log_start)).exp();
            times.push(value);
        }
        times.push(t_end);

        times
    }
}


