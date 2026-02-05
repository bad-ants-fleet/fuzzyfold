use clap::Parser;
use rand::rng;
use rand::prelude::IndexedRandom;

/// Generate random sequences from a given alphabet.
#[derive(Parser, Debug)]
#[command(author, version, about = "Generate random sequences", long_about = None)]
struct Args {
    /// Alphabet to choose from, comma-separated (e.g., A,C,G,U)
    #[arg(short, long, default_value = "A,C,G,U")]
    alphabet: String,

    /// Length of each generated sequence
    #[arg(short, long, default_value_t = 50)]
    length: usize,

    /// Number of sequences to generate
    #[arg(short, long, default_value_t = 1)]
    num: usize,

    /// Print the open chain configuration under each sequence
    #[arg(short, long)]
    eval: bool,

}

fn main() {
    let args = Args::parse();

    // Parse alphabet
    let alphabet: Vec<char> = args.alphabet.split(',').map(|s| s.chars().next().unwrap()).collect();

    // Prepare RNG
    let mut rng = rng();

    for _ in 0..args.num {
        let seq: String = (0..args.length)
            .map(|_| *alphabet.choose(&mut rng).unwrap())
            .collect();
        println!("{}", seq);
        if args.eval {
            println!("{}", ".".repeat(args.length));
        }
    }
}
