# `findpath_pseudo` — Performance Bottleneck Analysis

Profiled 2026-08-04 using `samply` + Criterion on macOS (Apple Silicon).
Workload: synthetic H-type PK structures, 30–90 nt, beam widths 1–200.

## Criterion benchmark results (release build)

### Beam-width sweep — 50 nt / 20 pairs

| beam | before  | after   | change |
|------|---------|---------|--------|
| 1    | 992 µs  | 604 µs  | −39%   |
| 5    | 4.4 ms  | 2.67 ms | −40%   |
| 10   | 8.3 ms  | 5.13 ms | −38%   |
| 25   | 19 ms   | 11.6 ms | −39%   |
| 50   | 36 ms   | 22.2 ms | −38%   |

Scaling remains **linear in beam width** post-optimisation.

### Structure-size sweep — beam = 10

| structure         | before | after   | change |
|-------------------|--------|---------|--------|
| 30 nt / 10 pairs  | 1.1 ms | 669 µs  | −42%   |
| 50 nt / 20 pairs  | 8.3 ms | 5.13 ms | −38%   |
| 70 nt / 30 pairs  | 28 ms  | 17.4 ms | −38%   |
| 90 nt / 40 pairs  | 67 ms  | 40.4 ms | −40%   |

### Non-empty start vs. empty start — 70 nt, beam = 10

| start            | before  | after   | change |
|------------------|---------|---------|--------|
| empty            | 28.5 ms | 17.4 ms | −39%   |
| stem-2 preformed | 8.9 ms  | 5.70 ms | −36%   |

The ~39% gain is uniform across structure sizes and beam widths: B1
(eliminating the O(P²) `pair_table_to_dot_bracket` + full `parse_structure`
re-parse on every cache miss) is the dominant contributor.  B2–B4 account
for the remainder and are especially visible at large beam widths.

---

## Bottleneck 1 — PairTable → String → PairTable roundtrip (highest priority)

**Location:** `pseudo_findpath.rs::eval_energy` calling into
`ff_energy::pair_table_to_dot_bracket` + `ff_energy::parse_structure`.

Every cache-miss energy evaluation does:

```
PairTable
  → pair_table_to_dot_bracket    (neighbors.rs:91)   O(P²) crossing-graph coloring
  → String
  → parse_structure              (enumerate.rs:225)
      → ExtendedDotBracketVec::try_from              O(N) lex
      → extended_dot_bracket_to_pair_table           O(N) + new HashMap + new PairTable
      → build_closed_regions_tree                    O(N) amortized
      → enumerate_loops                              O(N) + O(B·L) for PK bands
```

`pair_table_to_dot_bracket` is O(P²) (inner loop at `neighbors.rs:110`).
For 40 pairs it does ~800 comparisons just to produce a string that
`parse_structure` immediately re-parses back into a new PairTable.
The original PairTable was already in hand.

**Fix:** add `parse_loops_from_pt(pt: &PairTable) -> Vec<Loop>` in `ff_energy`
that calls `build_closed_regions_tree(pt)` + `enumerate_loops` directly,
skipping the string entirely.  `pair_table_to_dot_bracket` then only runs
once per winning step for human-readable output, not per candidate.

**Expected gain:** ~2–3× per unique eval (eliminates the largest fraction of
work in the hot path).

---

## Bottleneck 2 — Premature `path` / `remaining_moves` cloning (highest priority)

**Location:** `pseudo_findpath.rs:214–218`, inside the inner expansion loop.

```rust
let mut new_remaining = parent.remaining_moves.clone();  // O(moves)
new_remaining.remove(idx);
let mut new_path = parent.path.clone();                  // O(path_len)
new_path.push(mv.clone());
```

These clones run for **every candidate** before dedup and beam truncation.
At step k with beam=50 and 40 total moves:
- Up to `50 × (40 - k) ≈ 1000` candidates generated per step
- After dedup + truncation: only 50 survive
- ~950 of the 1000 `path` + `remaining_moves` clones are wasted

**Fix:** restructure the expansion loop so it only stores lightweight values
`(new_pt, new_saddle, energy, parent_index, move_index)` during expansion.
After dedup + truncation, reconstruct the full `remaining_moves` and `path`
for the `beam_width` survivors only.

**Expected gain:** ~10–20× fewer clone operations per step; dominant effect
at large beam widths.

---

## Bottleneck 3 — SipHash for the energy cache (easy win)

**Location:** `pseudo_findpath.rs:135`, `HashMap<PairTable, f64>`.

Rust's default SipHash hashes 180 bytes (90 nt × 2 bytes/entry) per lookup.
The cache is probed `beam × remaining` times per step — ~80 K probes for
beam=50, 90 nt.  `FxHashMap` from `rustc-hash` (already in `Cargo.toml` as
a workspace dependency) is ~3× faster for non-adversarial integer-heavy keys.

**Fix:** one import + one type annotation change:

```rust
use rustc_hash::FxHashMap;
let mut cache: FxHashMap<PairTable, f64> = FxHashMap::default();
```

**Expected gain:** ~3× faster per cache probe; meaningful at large beam widths.

---

## Bottleneck 4 — `Vec::remove` is O(moves) (easy win)

**Location:** `pseudo_findpath.rs:215`.

```rust
new_remaining.remove(idx);  // shifts all elements after idx — O(moves)
```

Move ordering within `remaining_moves` doesn't affect correctness (the beam
re-iterates the whole list fresh at every step).  `swap_remove(idx)` replaces
the removed element with the last element in O(1).

**Fix:** `new_remaining.swap_remove(idx);`

**Expected gain:** O(1) instead of O(N) per removal; low absolute impact
but trivially free.

---

## Non-issue: `build_closed_regions_tree`

The tree build (`closed_region_tree.rs:72`) scans N positions once (O(N))
and calls `add_to_tree` which sorts children by `i`.  For H-type PKs the
region count is small (≤2 top-level regions, 0 children), so the sort is
over at most 2 elements.  This is **not** a bottleneck — the cost is
dominated by the string roundtrip (Bottleneck 1) that feeds into it.

---

## Implementation order

| Priority | Change | File(s) | Effort |
|----------|--------|---------|--------|
| 1 | Add `parse_loops_from_pt(pt)` | `ff_energy/src/pseudoknots/enumerate.rs` + `mod.rs` | Medium |
| 1 | Use it in `eval_energy` | `ff_findpath/src/pseudo_findpath.rs` | Trivial |
| 2 | Defer `path`/`remaining_moves` cloning | `ff_findpath/src/pseudo_findpath.rs` | Medium |
| 3 | `FxHashMap` for cache | `ff_findpath/src/pseudo_findpath.rs` | Trivial |
| 4 | `swap_remove` | `ff_findpath/src/pseudo_findpath.rs` | Trivial |

Run `cargo bench -p ff_findpath --bench bench_pseudo_findpath` before and after
each change to track the effect in isolation.
