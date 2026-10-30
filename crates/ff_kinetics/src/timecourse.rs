//! Macrostate counts of many trajectories at common output times, with the
//! file formats .tln (counts, JSON) and .nxy (occupancies), and a quick plot.
//!
//! A trajectory is the macrostate index (of a MacrostateRegistry) at each
//! output time. The registry gives the names and lengths when reading,
//! writing or plotting.

use std::fmt;
use std::fs;
use std::error::Error;
use std::io::Write;
use std::path::Path;
use serde::{Serialize, Deserialize};
use plotters::prelude::*;
use plotters::style::Palette99;

use ff_energy::EnergyModel;
use crate::macrostates::MacrostateRegistry;

#[derive(Debug)]
pub enum TimecourseError {
    Io(std::io::Error),
    Json(serde_json::Error),
    TimepointCountMismatch { found: usize, expected: usize },
    TimeMismatch { file_time: f64, expected_time: f64 },
    MacrostateNotFound(String),
}

impl fmt::Display for TimecourseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {}", e),
            Self::Json(e) => write!(f, "JSON parse error: {}", e),
            Self::TimepointCountMismatch { found, expected } =>
                write!(f, "Timecourse file has {found} timepoints, expected {expected}"),
            Self::TimeMismatch { file_time, expected_time } =>
                write!(f, "Time mismatch: {file_time} vs {expected_time}"),
            Self::MacrostateNotFound(name) =>
                write!(f, "Macrostate '{name}' not found in registry"),
        }
    }
}

impl Error for TimecourseError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::Json(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for TimecourseError {
    fn from(e: std::io::Error) -> Self { Self::Io(e) }
}

impl From<serde_json::Error> for TimecourseError {
    fn from(e: serde_json::Error) -> Self { Self::Json(e) }
}

/// Counts of macrostates (by registry index) at each output time.
#[derive(Debug, Clone, PartialEq)]
pub struct Timecourse {
    pub times: Vec<f64>,
    /// `counts[t][m]`: the number of trajectories in macrostate m at `times[t]`.
    pub counts: Vec<Vec<usize>>,
}

/// The .tln file: per time point, (length, name, count) of the occupied
/// macrostates and the number of trajectories.
#[derive(Serialize, Deserialize)]
struct TlnFile {
    points: Vec<TlnPoint>,
}

#[derive(Serialize, Deserialize)]
struct TlnPoint {
    time: f64,
    ensemble: Vec<(usize, String, usize)>,
    counter: usize,
}

fn same_time(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

impl Timecourse {
    /// No trajectories at the given times, for num_macrostates macrostates.
    pub fn new(times: &[f64], num_macrostates: usize) -> Self {
        Self {
            times: times.to_vec(),
            counts: vec![vec![0; num_macrostates]; times.len()],
        }
    }

    /// Add a trajectory: its macrostate index at each output time.
    pub fn add_trajectory(&mut self, macrostates: &[usize]) {
        assert_eq!(macrostates.len(), self.times.len(), "one macrostate per output time");
        for (row, &m) in self.counts.iter_mut().zip(macrostates) {
            row[m] += 1;
        }
    }

    /// Add the counts of another time course with the same times.
    pub fn merge(&mut self, other: &Timecourse) {
        assert_eq!(self.times, other.times, "cannot merge different output times");
        for (row, other_row) in self.counts.iter_mut().zip(&other.counts) {
            assert_eq!(row.len(), other_row.len(), "cannot merge different macrostates");
            for (c, o) in row.iter_mut().zip(other_row) {
                *c += o;
            }
        }
    }

    /// The number of trajectories at time index t.
    pub fn num_trajectories(&self, t: usize) -> usize {
        self.counts[t].iter().sum()
    }

    /// Macrostate names in registry order, and the counts per name (summed
    /// over the lengths of co-transcriptional macrostates): `[t][name]`.
    pub fn counts_by_name<'a, E: EnergyModel>(
        &self,
        registry: &'a MacrostateRegistry<E>,
    ) -> (Vec<&'a str>, Vec<Vec<usize>>) {
        let mut names: Vec<&str> = Vec::new();
        let mut column: Vec<usize> = Vec::new();
        for (_, (_, ms)) in registry.iter() {
            let i = names.iter().position(|n| *n == ms.name()).unwrap_or_else(|| {
                names.push(ms.name());
                names.len() - 1
            });
            column.push(i);
        }
        let counts = self.counts.iter().map(|row| {
            let mut by_name = vec![0; names.len()];
            for (m, &c) in row.iter().enumerate() {
                by_name[column[m]] += c;
            }
            by_name
        }).collect();
        (names, counts)
    }

    /// Read a .tln file, which must have the given output times and only
    /// macrostates of the registry.
    pub fn read_tln<E: EnergyModel>(
        path: impl AsRef<Path>,
        times: &[f64],
        registry: &MacrostateRegistry<E>,
    ) -> Result<Self, TimecourseError> {
        let file: TlnFile = serde_json::from_str(&fs::read_to_string(path)?)?;
        if file.points.len() != times.len() {
            return Err(TimecourseError::TimepointCountMismatch {
                found: file.points.len(),
                expected: times.len(),
            });
        }
        let mut tc = Timecourse::new(times, registry.len());
        for (row, (point, &time)) in tc.counts.iter_mut().zip(file.points.iter().zip(times)) {
            if !same_time(point.time, time) {
                return Err(TimecourseError::TimeMismatch {
                    file_time: point.time,
                    expected_time: time,
                });
            }
            for (length, name, count) in &point.ensemble {
                let (m, _) = registry.iter()
                    .find(|(_, (l, ms))| l == length && ms.name() == name)
                    .ok_or_else(|| TimecourseError::MacrostateNotFound(name.clone()))?;
                row[m] += count;
            }
        }
        Ok(tc)
    }

    /// Write a .tln file.
    pub fn write_tln<E: EnergyModel>(
        &self,
        path: impl AsRef<Path>,
        registry: &MacrostateRegistry<E>,
    ) -> Result<(), TimecourseError> {
        let macrostates = registry.macrostates();
        let points = self.times.iter().zip(&self.counts).map(|(&time, row)| TlnPoint {
            time,
            ensemble: row.iter().enumerate()
                .filter(|(_, c)| **c > 0)
                .map(|(m, &c)| (macrostates[m].0, macrostates[m].1.name().to_string(), c))
                .collect(),
            counter: row.iter().sum(),
        }).collect();
        fs::write(path, serde_json::to_string_pretty(&TlnFile { points })?)?;
        Ok(())
    }

    /// Write a .nxy file: the occupancy of each macrostate name over time.
    pub fn write_nxy<E: EnergyModel>(
        &self,
        writer: &mut impl Write,
        registry: &MacrostateRegistry<E>,
    ) -> std::io::Result<()> {
        let (names, counts) = self.counts_by_name(registry);
        write!(writer, "{:>13}", "time")?;
        for name in &names {
            write!(writer, " {:>13}", name)?;
        }
        writeln!(writer)?;
        for (time, row) in self.times.iter().zip(&counts) {
            let total = row.iter().sum::<usize>().max(1) as f64;
            write!(writer, "{:13.6e}", time)?;
            for &c in row {
                write!(writer, " {:13.6e}", c as f64 / total)?;
            }
            writeln!(writer)?;
        }
        Ok(())
    }

    /// A quick SVG plot: linear up to t_sep, logarithmic up to t_end, with
    /// the standard errors as bands. Shows the first macrostate
    /// (Unassigned) and those with an occupancy of at least 0.02.
    pub fn plot_svg<E: EnergyModel>(
        &self,
        path: impl AsRef<Path>,
        registry: &MacrostateRegistry<E>,
        title: &str,
        t_sep: f64,
        t_end: f64,
    ) -> Result<(), Box<dyn Error + Send + Sync>> {
        if !(t_sep > 0.0 && t_end > t_sep) {
            return Err(format!("plot needs 0 < t_sep ({t_sep}) < t_end ({t_end})").into());
        }

        let root = SVGBackend::new(path.as_ref(), (1024, 480)).into_drawing_area();
        root.fill(&WHITE)?;
        root.titled(title, ("sans-serif", 28))?;
        root.draw_text(
            "time",
            &("sans-serif", 22).into_font().into_text_style(&root),
            (496, 450),
        )?;

        let eps = 1e-12;
        let (left, right) = root.split_horizontally(512);

        let mut chart_left = ChartBuilder::on(&left)
            .caption("Linear plot", ("sans-serif", 18))
            .margin(20)
            .margin_top(40)
            .margin_right(0)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(0.0..(t_sep + eps), 0.0..1.0)?;
        chart_left
            .configure_mesh()
            .y_desc("occupancy")
            .x_labels(6)
            .x_label_formatter(&|x| if *x < 0.01 { format!("{:.1e}", x) } else { format!("{}", x) })
            .y_labels(10)
            .light_line_style(TRANSPARENT)
            .axis_desc_style(("sans-serif", 22))
            .label_style(("sans-serif", 18))
            .draw()
            ?;
        chart_left.draw_series(std::iter::once(PathElement::new(
            vec![(t_sep, 0.0), (t_sep, 1.0)],
            BLACK.mix(0.7),
        )))?;

        let mut chart_right = ChartBuilder::on(&right)
            .caption("Logarithmic plot", ("sans-serif", 18))
            .margin(20)
            .margin_top(40)
            .margin_left(0)
            .margin_right(40)
            .x_label_area_size(40)
            .y_label_area_size(0)
            .build_cartesian_2d(((t_sep + eps)..(t_end + eps)).log_scale(), 0.0..1.0)
            ?;
        chart_right
            .configure_mesh()
            .x_labels(6)
            .x_label_formatter(&|x| if *x < 0.01 { format!("{:.1e}", x) } else { format!("{}", x) })
            .y_labels(10)
            .light_line_style(TRANSPARENT)
            .label_style(("sans-serif", 18))
            .draw()?;
        chart_right.draw_series(std::iter::once(PathElement::new(
            vec![(t_sep, 0.0), (t_sep, 1.0)],
            BLACK.mix(0.7),
        )))?;

        // (time, occupancy, standard error) per name.
        let (names, counts) = self.counts_by_name(registry);
        let mut shown = Vec::new();
        for (i, name) in names.iter().enumerate() {
            let series: Vec<(f64, f64, f64)> = self.times.iter().zip(&counts).map(|(&t, row)| {
                let n = row.iter().sum::<usize>();
                if n == 0 {
                    return (t, 0.0, 0.0);
                }
                let p = row[i] as f64 / n as f64;
                (t, p, (p * (1.0 - p) / n as f64).sqrt())
            }).collect();
            if i == 0 || series.iter().any(|(_, p, _)| *p >= 0.02) {
                shown.push((name, series));
            }
        }

        for (i, (name, series)) in shown.iter().enumerate() {
            let color = Palette99::pick(i).mix(0.9);
            let band_color = color.mix(0.2);
            let upper = series.iter().map(|(t, p, se)| (*t, (p + se).min(1.0)));
            let lower = series.iter().rev().map(|(t, p, se)| (*t, (p - se).max(0.0)));
            let band = upper.chain(lower);

            chart_left.draw_series(AreaSeries::new(
                band.clone().filter(|(t, _)| *t <= t_sep),
                0.0,
                band_color,
            ))?;
            chart_left.draw_series(LineSeries::new(
                series.iter().map(|(t, p, _)| (*t, *p)).filter(|(t, _)| *t <= t_sep + eps),
                color.stroke_width(2),
            ))?;
            chart_right.draw_series(AreaSeries::new(
                band.filter(|(t, _)| *t >= t_sep),
                0.0,
                band_color,
            ))?;
            chart_right.draw_series(LineSeries::new(
                series.iter().map(|(t, p, _)| (*t, *p)).filter(|(t, _)| *t >= t_sep - eps),
                color.stroke_width(2),
            ))?
                .label(name.to_string())
                .legend(move |(x, y)| {
                    PathElement::new(vec![(x, y), (x + 20, y)], color.stroke_width(2))
                });
        }

        chart_right
            .configure_series_labels()
            .border_style(BLACK)
            .background_style(WHITE.mix(0.8))
            .position(SeriesLabelPosition::UpperRight)
            .label_font(("sans-serif", 16).into_font())
            .draw()?;
        root.present()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;
    use std::io::Cursor;
    use ff_energy::{NucleotideVec, Fuzzyfold};

    fn registry() -> MacrostateRegistry<Fuzzyfold> {
        let seq = NucleotideVec::try_from("GGGAAACCC").unwrap();
        let mut registry = MacrostateRegistry::from((seq, Arc::new(Fuzzyfold::default())));
        let ms = b">open\nGGGAAACCC\n.........\n";
        registry.insert_from_reader(Cursor::new(ms), "test", 0).unwrap();
        let ms = b">hairpin\nGGGAAACCC\n(((...)))\n";
        registry.insert_from_reader(Cursor::new(ms), "test", 0).unwrap();
        registry
    }

    #[test]
    fn counts_merge_and_tln_round_trip() {
        let registry = registry();
        let times = [0.0, 0.5, 1.0];
        let mut a = Timecourse::new(&times, registry.len());
        a.add_trajectory(&[1, 1, 2]);
        a.add_trajectory(&[1, 0, 2]);
        let mut b = Timecourse::new(&times, registry.len());
        b.add_trajectory(&[1, 2, 2]);
        a.merge(&b);
        assert_eq!(a.counts, vec![vec![0, 3, 0], vec![1, 1, 1], vec![0, 0, 3]]);
        assert_eq!(a.num_trajectories(1), 3);

        let dir = std::env::temp_dir().join(format!("ff_tc_{}.tln", std::process::id()));
        a.write_tln(&dir, &registry).unwrap();
        let read = Timecourse::read_tln(&dir, &times, &registry).unwrap();
        fs::remove_file(&dir).unwrap();
        assert_eq!(read, a);
        assert!(matches!(Timecourse::read_tln("does/not/exist.tln", &times, &registry),
            Err(TimecourseError::Io(_))));

        let mut nxy = Vec::new();
        a.write_nxy(&mut nxy, &registry).unwrap();
        let nxy = String::from_utf8(nxy).unwrap();
        let lines: Vec<&str> = nxy.lines().collect();
        assert_eq!(lines[0].split_whitespace().collect::<Vec<_>>(),
            ["time", "Unassigned", "open", "hairpin"]);
        assert_eq!(lines[1].split_whitespace().collect::<Vec<_>>(),
            ["0.000000e0", "0.000000e0", "1.000000e0", "0.000000e0"]);
    }
}
