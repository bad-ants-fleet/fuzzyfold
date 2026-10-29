use ff_findpath::findpath::{findpath};
use ff_energy::ViennaRNA;

use std::fs;
use std::path::Path;
use std::io;

use clap::Parser;



/// Simple CLI for getting the input Path
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of file containing input data
    #[arg(short, long, value_name = "FILENAME",)]
    filename: String,

    /// Search width parameter for findpath (default= 1: corresponds to Greedy algorithm)
    #[arg(short = 'm', long = "search-width", value_name = "SEARCH-WIDTH", default_value = "1")]
    m: usize,

}

/// Reads sequence and structures from a three-line test file.
/// Returns a tuple (sequence, s1, s2) on success.
pub fn load_test_file(path: &str) -> Result<(String, String, String), io::Error> {
    let file_path = Path::new(path);

    // 1. Check if file exists (matches `if not file_path.is_file()`)
    if !file_path.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("Test file '{}' does not exist", path),
        ));
    }

    // 2. Read file content
    let content = fs::read_to_string(file_path)?;

    // 3. Process lines: strip whitespace and remove empty lines
    let lines: Vec<&str> = content
        .lines()
        .map(|line| line.trim())
        .filter(|line| !line.is_empty())
        .collect();

    // 4. Validate line count
    if lines.len() < 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Expected at least three non-empty lines in {}", path),
        ));
    }

    // 5. Return &str slices for sequence and structures
    Ok((
        lines[0].to_string(),
        lines[1].to_string(),
        lines[2].to_string(),
    ))
}


fn main()  -> Result<(), Box<dyn std::error::Error>> {
    // Example usage of load_test_file

    let args = Args::parse();
    let test_file_path = args.filename;

    // call the sequence, structures from the test file
    let (sequence, struct1, struct2) = load_test_file(&test_file_path)?; // These are Strings
    let seq: &str = &sequence; // Convert String to &str
    let struct1: &str = &struct1; // Convert String to &str
    let struct2: &str = &struct2;
    //let pt1 = PairTable::try_from(struct1).unwrap(); 
    //let pt2 = PairTable::try_from(struct2).unwrap();

    // Initialize Model (This might fail if params aren't found in your env)
    // Ensure you have valid params or a mock model available.
    let model = ViennaRNA::default();

    let (steps, stats) = findpath(&model, seq, struct1, struct2, args.m, None).unwrap();
    println!("Folding Path found with findpath algorithm:");
    println!("-----------------");
    println!("{} \t applied move \t energy", seq);
    for step in steps {
        println!("{} \t {} \t {} kcal/mol", step.structure, step.move_applied.unwrap_or_default(), step.energy );
    }   
    println!("-----------------");
    println!("Statistics:");
    println!("{}", stats);
    Ok(())

}