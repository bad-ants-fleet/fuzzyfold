use domainlevel::Acfp;
use domainlevel::DomainRegistry;
use domainlevel::SegmentSequence;

fn main() {
    let mut registry = DomainRegistry::new();

    //let acfp = Acfp::try_from(". () (). .(.) .(..) .((.))").unwrap();
    //let acfp = Acfp::try_from(". .. .() ..() (.()) ((()))").unwrap();
    //let acfp = Acfp::try_from(". () (). .(.) .(.). .((.))").unwrap();
    let acfp = Acfp::try_from(". () (). .(.) .(..) .(...) .(...). .((...))").unwrap();

    println!("Is cycle free: {}", acfp.is_cycle_free());
    println!("Has valid pair hierarchy: {}", acfp.has_valid_pair_hierarchy());
    println!("Pair hierarchy: {:?}", acfp.pair_hierarchy());
    println!("Has valid sequence design: {}", acfp.has_valid_design(&mut registry));

    println!();
    let segseq = SegmentSequence::design_from_acfp(&acfp, &mut registry).unwrap();
    println!("{}", &segseq.get_domain_sequence().iter()
        .map(|d| format!("{}", d)).collect::<Vec<_>>().join(" "));

    //for db in &segseq.get_domain_transcription_complexes(&registry) {
    //    println!("{:<2} {}", db.len(), db);
    //}
    for db in &segseq.get_segment_transciption_complexes(&registry) {
        println!("{:<2} {}", db.len(), db);
    }
    for db in &segseq.get_adibatic_acfp(&registry) {
        println!("{:<2} {}", db.len(), db);
    }

    assert!(segseq.implements_acfp(acfp.path(), &registry));
}


