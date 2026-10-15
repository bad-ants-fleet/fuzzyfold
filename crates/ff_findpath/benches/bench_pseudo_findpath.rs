use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use std::hint::black_box;

use ff_energy::{parameters::{RNA_DP09, RNA_MT09}, ViennaRNA};
use ff_findpath::pseudo_findpath::findpath_pseudo;

// ---------------------------------------------------------------------------
// Test structures — synthetic H-type pseudoknots of increasing size.
//
// For n_stem pairs each helix, total length = 4·n_stem + loop_len.
// Loop length fixed at 10 throughout.
//
// Structure layout (n_stem = N):
//   pos  0 ..  N-1 : stem-1 opens  "((((("
//   pos  N .. 2N-1 : stem-2 opens  "[[[[["
//   pos 2N .. 3N-1 : stem-1 closes "))))))"
//   pos 3N .. 3N+9 : loop          ".........."
//   pos 3N+10 .. 4N+9 : stem-2 closes "]]]]]]"
// ---------------------------------------------------------------------------

struct Case {
    label:  &'static str,
    seq:    &'static str,
    target: &'static str,
}

static CASES: &[Case] = &[
    Case {
        label:  "30nt/10pairs",
        seq:    "GGGGGGGGGGCCCCCAAAAAAAAAACCCCC",
        //       (((((  [[[[[  )))))  loop10  ]]]]]
        target: "((((([[[[[)))))..........]]]]]",
    },
    Case {
        label:  "50nt/20pairs",
        seq:    "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCAAAAAAAAAACCCCCCCCCC",
        target: "(((((((((([[[[[[[[[[))))))))))..........]]]]]]]]]]",
    },
    Case {
        label:  "70nt/30pairs",
        seq:    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCC",
        target: "((((((((((((((([[[[[[[[[[[[[[[)))))))))))))))..........]]]]]]]]]]]]]]]",
    },
    Case {
        label:  "90nt/40pairs",
        seq:    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCAAAAAAAAAACCCCCCCCCCCCCCCCCCCC",
        target: "(((((((((((((((((((([[[[[[[[[[[[[[[[[[[[))))))))))))))))))))..........]]]]]]]]]]]]]]]]]]]]",
    },
];

fn make_model() -> ViennaRNA {
    ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09)
}

// ---------------------------------------------------------------------------
// 1. Beam-width sweep on a fixed structure (50 nt)
// ---------------------------------------------------------------------------
fn bench_beam_width(c: &mut Criterion) {
    let model = make_model();
    let case = &CASES[1]; // 50 nt / 20 pairs

    let mut group = c.benchmark_group("pseudo_findpath/beam_width");
    for &beam in &[1usize, 5, 10, 25, 50] {
        group.bench_with_input(BenchmarkId::from_parameter(beam), &beam, |b, &beam| {
            b.iter(|| {
                findpath_pseudo(
                    black_box(&model),
                    black_box(case.seq),
                    black_box(None),
                    black_box(case.target),
                    black_box(beam),
                    black_box(None),
                    false,
                )
                .unwrap()
            })
        });
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 2. Structure-size sweep at fixed beam width (10)
// ---------------------------------------------------------------------------
fn bench_structure_size(c: &mut Criterion) {
    let model = make_model();
    const BEAM: usize = 10;

    let mut group = c.benchmark_group("pseudo_findpath/structure_size");
    for case in CASES {
        group.bench_with_input(
            BenchmarkId::from_parameter(case.label),
            case,
            |b, case| {
                b.iter(|| {
                    findpath_pseudo(
                        black_box(&model),
                        black_box(case.seq),
                        black_box(None),
                        black_box(case.target),
                        black_box(BEAM),
                        black_box(None),
                        false,
                    )
                    .unwrap()
                })
            },
        );
    }
    group.finish();
}

// ---------------------------------------------------------------------------
// 3. Non-empty start vs empty start (70 nt, beam=10)
// ---------------------------------------------------------------------------
fn bench_start_effect(c: &mut Criterion) {
    let model = make_model();
    let case = &CASES[2]; // 70 nt / 30 pairs

    // Start = second stem pre-formed (only 15 insertions for stem-1 needed).
    // 70 nt: stem-2 left at pos 15-29, right at pos 55-69.
    // Constructed as: 15 dots + 15 opens + 25 dots + 15 closes = 70 chars.
    // 15 dots + 15 '(' + 25 dots + 15 ')' = 70 chars
    let start_preformed = concat!(
        "...............",   // 15 dots  (stem-1 left, unpaired in start)
        "(((((((((((((((", // 15 opens (stem-2 left)
        ".........................", // 25 dots (stem-1 right + loop, unpaired in start)
        ")))))))))))))))");  // 15 closes (stem-2 right)

    let mut group = c.benchmark_group("pseudo_findpath/start_effect");
    group.bench_function("empty_start", |b| {
        b.iter(|| {
            findpath_pseudo(
                black_box(&model),
                black_box(case.seq),
                black_box(None),
                black_box(case.target),
                black_box(10usize),
                black_box(None),
                false,
            )
            .unwrap()
        })
    });
    // Only run if lengths match (compile-time check via assert at top-level would
    // turn a typo into a test failure rather than a silently-skipped bench).
    assert_eq!(start_preformed.len(), case.seq.len(),
        "start_preformed length mismatch for bench_start_effect");
    group.bench_function("preformed_stem2", |b| {
        b.iter(|| {
            findpath_pseudo(
                black_box(&model),
                black_box(case.seq),
                black_box(Some(start_preformed)),
                black_box(case.target),
                black_box(10usize),
                black_box(None),
                false,
            )
            .unwrap()
        })
    });
    group.finish();
}

criterion_group!(benches, bench_beam_width, bench_structure_size, bench_start_effect);
criterion_main!(benches);
