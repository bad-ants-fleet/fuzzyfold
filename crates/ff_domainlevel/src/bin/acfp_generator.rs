use itertools::Itertools;
use structure::DotBracketVec;
use domainlevel::Acfp;
use domainlevel::DomainRegistry;
use domainlevel::SegmentSequence;
use domainlevel::{enum_structs, generate_structs};

fn display_acfp(acfp: &[DotBracketVec]) -> String {
    acfp.iter()
        .map(|s| s.iter()
            .map(|&db| char::from(db)) 
            .collect::<String>())
        .join(" ")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {

    for n in 1..=10 {
        let (val, _) = enum_structs(n);
        println!("Length {n:>2}: {val} strucures.");
    }
    println!("--");

    for n in 1..=10 {
        let (_, seq) = enum_structs(n);
        let prod: usize = seq.iter().product();
        println!("Length {n:>2}: {prod} total paths.");
    }
    println!("--");

    let mut registry = DomainRegistry::new();

    let mut acfps: Vec<Acfp> = vec![Acfp::try_from(".").expect("must work")];
    for len in 2..=10 {
        let mut next_acfps = Vec::new();
        let struct_set = generate_structs(len);

        let mut lin_ext = 0;
        for acfp in &acfps {
            for db in &struct_set {
                let mut new = acfp.clone();
                new.extend_by_one(DotBracketVec(db.clone()));
                if new.is_valid(&mut registry) {
                    next_acfps.push(new.clone());
                    lin_ext += new.all_total_orders().unwrap().len();
                    if len == 6 {
                        let segseq = SegmentSequence::design_from_acfp(&new, &mut registry).unwrap();
                        println!("{len:<3} \"{}\" {}", 
                            display_acfp(new.path()), 
                            &segseq.get_domain_sequence().iter()
                                .map(|d| format!("{}", d)).collect::<Vec<_>>().join(" ")
                        );
                    }
                }
            }
        }
        println!("Length {:>2}: {} valid paths, {} linear extensions.", len, next_acfps.len(), lin_ext);
        acfps = next_acfps;
    }
    println!("--");

    Ok(())
}

