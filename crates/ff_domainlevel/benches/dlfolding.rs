use criterion::Criterion;
use criterion::criterion_group;
use criterion::criterion_main;

use ff_domainlevel::DomainRegistry;
use ff_domainlevel::NussinovDP;

pub fn dl_nussinov(c: &mut Criterion) {
    let mut group = c.benchmark_group("DomainLevel");
    //group.measurement_time(std::time::Duration::from_secs(50)); // increase from default 5s

    let mut registry = DomainRegistry::new();
    registry.intern("a", 1);
    let ndp = NussinovDP::try_from(("a a* a a a* a a* a a* a a*", &registry)).unwrap();

    group.bench_function("Enumerate all MFE structures.", |b| {
        b.iter(|| {
            let _ = ndp.all_mfe_structs(None);
        });
    });
}

criterion_group!(benches, dl_nussinov);
criterion_main!(benches);

