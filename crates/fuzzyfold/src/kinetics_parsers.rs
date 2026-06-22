use clap::Args;
use anyhow::bail;
use anyhow::Result;
use ff_energy::NucleotideVec;
use ff_kinetics::Arrhenius;
use ff_structure::DotBracketVec;
use plotters::prelude::*;

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


#[derive(Debug, Args)]
pub struct TimelineParameters {
    /// Extension time during transcription.
    #[arg(long)]
    pub t_ext: Option<f64>,   

    /// Simulation stop time/ Posttranscriptional simulation time (in cotranscriptional simulations).
    #[arg(long, default_value_t = 1.0)]
    pub t_end: f64,

    /// Number of time points on the linear scale.
    #[arg(long, default_value_t = 100)]
    pub t_lin: usize,

    /// Number of time points on the logarithmic scale.
    #[arg(long, default_value_t = 100)]
    pub t_log: usize,

    /// Which mode? t_sep given: seperator between linear and logarithimic part at t_sep, 
    /// t-lin: timepoints recorded on a linear timescale between 0 and t-sep, t-log: timepoints 
    /// recorded on a logarithmic timescale between t-sep and total time.
    #[arg(long)]
    pub t_sep: Option<f64>,   
}

impl TimelineParameters {
    /// Validate that all parameters make sense.
    pub fn validate(&self) -> Result<()> {
        
        if self.t_ext.is_none() { 

            if let Some(sep) = self.t_sep {

                if self.t_end <= sep {
                    bail!("t_end ({}) must be greater than t_sep ({})", self.t_end, sep);
                }
            }
        }
        
        if self.t_lin == 0 && self.t_log > 1 {
            bail!("t_lin must be > 0 if t_log > 1 (got t_lin={}, t_log={})", self.t_lin, self.t_log);
        }
        Ok(())
    }

    pub fn get_output_times(&self, sequence: &NucleotideVec, structure: Option<&DotBracketVec>) -> Vec<f64> {
        let t_end = self.t_end;
        let t_ext = self.t_ext;
        let t_lin = self.t_lin;
        let t_log = self.t_log;
        let t_sep = self.t_sep;
        let mut times = vec![0.0];

        if t_ext.is_none() { // full length simulation
            
            // t_sep given => t_lin evenly spaced timepoints on a linear timescale between 0 and t-sep and 
            // t_log timepoints on a logarithmic timescale 

                // Linear segments: append `t_lin` evenly spaced points between 0...t_sep
                let start = *times.last().unwrap();
                let step = t_sep.unwrap() / t_lin as f64;
                for i in 1..=t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic tail: append 't_log logarithmic timepoints between t-sep...t_end
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

        } else {  // Co-transcriptional simulation 
            
            if t_sep.is_none() { 
                // no t_sep given => t_lin = timepoints per transcript length; t_log = timepoints selected for full length simulation 
                // Co-transcriptional folding: linear timescale, posttranscriptional folding: logarithmic timescale 

                // Linear time points: append `t_lin` evenly spaced points per t_ext  (0...t-ext)
                let start_len = structure.map(|s| s.len()).unwrap_or(1);
                let total_len = sequence.len();

                let mut len = start_len;

                while len < total_len {
                    let start = *times.last().unwrap();
                    let step = t_ext.unwrap() / t_lin as f64;
                    for i in 1..= t_lin {
                        times.push(start + i as f64 * step);
                    }
                    len += 1;
                }

                // Logarithmic tail: append 't-log' logarithmic spaced timepoints between t-ext * (sequence_len - start_len -1)...t_end
                let sep = *times.last().unwrap();
                let log_start = sep.ln();
                let log_end = (sep + t_end).ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(sep + t_end);

                times

            } else { 
            // t_sep given => t_lin evenly spaced timepoints on a linear timescale between 0 and t-sep and 
            // t_log timepoints on a logarithmic timescale between t-sep and (t-ext * (sequence_len - start_len -1) + t_end)

                let start_len = structure.map(|s| s.len()).unwrap_or(1);
                let total_len = sequence.len(); 

                // Linear segments: append `t_lin` evenly spaced points
                let start = *times.last().unwrap();
                let step = t_sep.unwrap() / t_lin as f64;
                for i in 1..= t_lin {
                    times.push(start + i as f64 * step);
                }

                // Logarithmic tail
                let start = *times.last().unwrap();
                let log_start = start.ln();
                let end = (t_ext.unwrap() * (total_len - start_len).as_f64()) + t_end;
                let log_end = end.ln();
                for i in 1..t_log {
                    let frac = i as f64 / t_log as f64;
                    let value = (log_start + frac * (log_end - log_start)).exp();
                    times.push(value);
                }
                times.push(end);

                times
            }
        }
    }
}


   
