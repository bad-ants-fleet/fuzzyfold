use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;

use ff_findpath::utils::{compare_structures, generate_valid_neighbors};
use ff_structure::PairTable;


fn bench_generate_valid_neighbors(c: &mut Criterion) {
    // ----- fixed benchmark inputs -----
    // Build a reasonable test PairTable and moves once, outside the iter loop.
    // This test is the large.txt from the test_data
    let current_pt = PairTable::try_from("......((...(((.(.(((((.......))))).))))....)).(((((((...))))....(((((....)))))..((.........(((((((((.(...((......))...).)))))).......))).........)))))").unwrap();
    let end_pt = PairTable::try_from("..(.......((((...(((((.......)))))))))(((.((...)).))).(((((...))))).)...((((((((.....))).)))))((((((.(...((......))...).)))))).......(((((.......)))))").unwrap();
 
    let diff = compare_structures(&current_pt, &end_pt);
    let available_moves = diff.move_list;

    
   

    c.bench_function("generate_valid_neighbors", |b| {
        b.iter(|| {
            // black_box to avoid over-optimization of inputs
            let pt = black_box(&current_pt);
            let moves = black_box(&available_moves[..]);

            // function under test
            let res = generate_valid_neighbors(pt, moves);

            // prevent compiler from discarding the result
            black_box(res);
        });
    });
}

criterion_group!(benches, bench_generate_valid_neighbors);
criterion_main!(benches);

