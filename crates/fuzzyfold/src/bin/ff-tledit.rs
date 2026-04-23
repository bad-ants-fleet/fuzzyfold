use colored::*;
use clap::Parser;
use anyhow::Result;
use std::sync::Arc;
use std::path::PathBuf;
use std::fs;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::io::BufWriter;
use serde_json::to_string_pretty;

use ff_energy::NucleotideVec;
use ff_kinetics::timeline::Timeline;
use ff_kinetics::timeline_plotting::plot_occupancy_over_time;
use ff_kinetics::MacrostateRegistry;

use fuzzyfold::energy_parsers::EnergyModelArguments;

#[derive(Debug, Parser)]
#[command(version, about = "Edit/import timeline format of simulated data.")]
pub struct Cli {
    #[arg(long, value_name = "FILE", num_args = 1.., required = false)]
    macrostates: Vec<PathBuf>,

    #[arg(long, num_args = 1..)]
    merge_tln: Vec<PathBuf>,

    #[arg(long, requires = "o_counter", num_args = 1..)]
    merge_nxy: Vec<PathBuf>,

    #[arg(long, requires = "merge_nxy")]
    o_counter: Option<usize>,

    #[arg(long, default_value_t = 1.0)]
    t_rescale: f64,

    #[arg(short, long)]
    title: String,

    #[arg(short, long, value_name = "FILE")]
    output: PathBuf,

    #[command(flatten, next_help_heading = "Energy model parameters")]
    energy: EnergyModelArguments,
}

fn get_sequence(msfile: &PathBuf) -> io::Result<NucleotideVec> {
    let file = File::open(msfile).expect("Failed to open first macrostate file");
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header_line = lines
        .next()
        .ok_or_else(|| io_err("Missing header line", &msfile.display().to_string()))??
        .trim()
        .to_string();
    let _ = header_line
        .strip_prefix('>')
        .ok_or_else(|| io_err("First line must start with '>'", &msfile.display().to_string()))?
        .trim()
        .to_string();

    let seq_line = lines
        .next()
        .ok_or_else(|| io_err("Missing sequence line", &msfile.display().to_string()))??;
    let sequence = NucleotideVec::try_from(seq_line.trim()).unwrap();
    Ok(sequence)
}

fn io_err(msg: &str, src: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("{} in {}", msg, src))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Build a Macrostate-registy from files
    let emodel = Arc::new(cli.energy.build_model());
    let sequence = get_sequence(&cli.macrostates[0]).unwrap();
    let mut registry = MacrostateRegistry::from((sequence.clone(), emodel));
    let _ = registry.insert_files(&cli.macrostates);
    // Verbose Output
    println!("# {:>4} {:<10} {} {:>5} {:>8}",
        "ID",
        "Macrostate".cyan(), format!("{}", sequence).yellow(), "Size", "Energy");
    for (id, m) in registry.iter() {
        if m.name() == "Unassigned" {
            println!("# {:4} {:<10}", 0, m.name());
            continue
        }
        println!("# {:4} {:<10} {:<} {:>5} {:>8.2}",
            id, 
            m.name(),
            m.get_lowest_microstate().unwrap(),
            m.len(),
            m.ensemble_energy().unwrap());
    }
    let shared_reg = Arc::new(registry);

    let tln_path = cli.output.with_extension("tln");
    let svg_path = cli.output.with_extension("svg");
    let nxy_path = cli.output.with_extension("nxy");

    // Build an empty timeline (NOTE: requires list of time points.)
    let mut timeline = Timeline::new(&[], Arc::clone(&shared_reg));

    // Merge existing timelines 
    for path in &cli.merge_nxy {
        timeline.load_nxy_data(path, cli.t_rescale, cli.o_counter
            .expect("needed for nxy"))?;
    }
    for path in &cli.merge_tln {
        timeline.load_tln_data(path, cli.t_rescale)?;
    }
    timeline.finalize();
    println!("{}", "Finished merging!".red());

    // save / print / plot.
    let mut writer = BufWriter::new(File::create(nxy_path.clone())?);
    write!(writer, "{}", timeline)?;
    println!("Wrote nxy file: {}", format!("{}",nxy_path.display()).green());
    plot_occupancy_over_time(&timeline, svg_path.clone(), &cli.title, 1e-5, 1e2);
    println!("Plotted svg file: {}", svg_path.display());
    let serial = timeline.to_serializable();
    let json = to_string_pretty(&serial).unwrap();
    fs::write(tln_path.clone(), json).unwrap();
    println!("Wrote tln file: {}", tln_path.display());

    Ok(())
}


