//! PK-aware findpath: beam search from an arbitrary start structure to a pseudoknotted target.
//!
//! Unlike the standard [`findpath`](crate::findpath::findpath), this variant:
//! - accepts any start structure (defaults to fully-unpaired if `None`)
//! - allows crossing base-pair insertions (pseudoknots)
//! - evaluates energy via [`energy_of_pseudoknotted_structure`] — a full
//!   closed-region-tree recalculation, memoized by pair table
//!
//! The algorithm is a directed beam search: only moves that close the diff
//! between start and target are considered (insertions of pairs missing in
//! start, deletions of pairs present in start but absent in target).

use std::cmp::Ordering;

use rustc_hash::{FxHashMap, FxHashSet};

use ff_energy::{
    pair_table_to_dot_bracket, parse_loops_from_pt, EnergyError, NucleotideVec,
    PseudoEnergyModel, ViennaRNA,
    extended_dot_bracket_to_pair_table,
};
use ff_structure::{ExtendedDotBracketVec, PairTable, StructureError};

use crate::utils::{Move, PathStats, PathStep, compare_structures};

// ---------------------------------------------------------------------------
// Internal beam state
// ---------------------------------------------------------------------------

struct PseudoIntermediate {
    pt:              PairTable,
    saddle_energy:   f64,
    remaining_moves: Vec<Move>,
    path:            Vec<Move>,
}

/// Lightweight record built during expansion — no `Vec` clones until survivors known (B2).
struct Expansion {
    pt:             PairTable,
    saddle_energy:  f64,
    current_energy: f64,
    parent_idx:     usize,
    move_idx:       usize,
}

// ---------------------------------------------------------------------------
// Energy helper (memoized)
// ---------------------------------------------------------------------------

/// Evaluate the full pseudoknot energy of `pt` in kcal/mol.
///
/// Results are stored in `cache` (keyed by pair table) so that the same
/// intermediate reached via different move orderings is only computed once.
/// Uses [`parse_loops_from_pt`] to skip the PT → string → PT round-trip (B1).
/// Returns `(energy, was_cache_hit)`.
fn eval_energy(
    model: &ViennaRNA,
    seq: &[ff_energy::Base],
    pt: &PairTable,
    cache: &mut FxHashMap<PairTable, f64>,
) -> Result<(f64, bool), String> {
    if let Some(&cached) = cache.get(pt) {
        return Ok((cached, true));
    }
    let loops = parse_loops_from_pt(pt);
    let energy = model.energy_of_pseudoknotted_structure(seq, &loops)
        .map(|e| e as f64 / 100.0)
        .map_err(|e: EnergyError| format!("energy: {e:?}"))?;
    cache.insert(pt.clone(), energy);
    Ok((energy, false))
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Find the minimum-barrier folding path from `start` to `target`.
///
/// # Arguments
/// - `model`       — ViennaRNA energy model (must have PK parameters attached via
///                   [`ViennaRNA::with_pseudoknot_params`] for PK-accurate energies)
/// - `sequence`    — RNA sequence string (`ACGU` alphabet)
/// - `start`       — start structure in extended dot-bracket notation, or `None`
///                   for the fully-unpaired structure
/// - `target`      — target structure in extended dot-bracket notation (e.g. `"((([[[)))...]]]"`)
/// - `beam_width`  — number of beam members kept per step (`1` = greedy)
/// - `max_energy`  — optional energy ceiling in kcal/mol; intermediates above
///                   this threshold are pruned
///
/// # Returns
/// `Ok((path, stats))` where `path` is the full step-by-step trajectory and
/// `stats` summarises the saddle energy and barrier.
///
/// # Note on bracket families
/// Structure strings in the returned path are produced by [`pair_table_to_dot_bracket`],
/// which applies a greedy graph-colouring. The family labels (`()`, `[]`, …) may differ
/// from those in `target` even when the pair tables are identical.
pub fn findpath_pseudo(
    model: &ViennaRNA,
    sequence: &str,
    start: Option<&str>,
    target: &str,
    beam_width: usize,
    max_energy: Option<f64>,
    cache_stats: bool,
) -> Result<(Vec<PathStep>, PathStats), String> {
    // ── parse inputs ─────────────────────────────────────────────────────────
    let seq_vec = NucleotideVec::try_from_rna(sequence)
        .map_err(|_| "Failed to parse RNA sequence".to_string())?;
    let seq: &[ff_energy::Base] = &seq_vec;

    let target_edbv = ExtendedDotBracketVec::try_from(target)
        .map_err(|e: StructureError| format!("Invalid target structure: {e}"))?;
    let target_pt = extended_dot_bracket_to_pair_table(&target_edbv)
        .map_err(|e: StructureError| format!("Target pair table: {e}"))?;

    let n = target_pt.len();
    if n != seq.len() {
        return Err(format!(
            "Sequence length ({}) != target structure length ({})", seq.len(), n
        ));
    }

    // ── start state ───────────────────────────────────────────────────────────
    let start_pt = match start {
        None => PairTable::new(n),
        Some(s) => {
            let edbv = ExtendedDotBracketVec::try_from(s)
                .map_err(|e: StructureError| format!("Invalid start structure: {e}"))?;
            let pt = extended_dot_bracket_to_pair_table(&edbv)
                .map_err(|e: StructureError| format!("Start pair table: {e}"))?;
            if pt.len() != n {
                return Err(format!(
                    "Start structure length ({}) != target structure length ({})", pt.len(), n
                ));
            }
            pt
        }
    };

    let start_dot = pair_table_to_dot_bracket(&start_pt)
        .unwrap_or_else(|_| ".".repeat(n));

    let mut cache: FxHashMap<PairTable, f64> = FxHashMap::default();
    let (start_energy, _) = eval_energy(model, seq, &start_pt, &mut cache)?;

    // ── move set ──────────────────────────────────────────────────────────────
    // compare_structures gives deletions-first, then insertions, sorted by i.
    // That ordering keeps deletions early so they unblock insertions in later steps.
    let diff = compare_structures(&start_pt, &target_pt);
    let total_steps = diff.move_list.len();

    if total_steps == 0 {
        let stats = PathStats {
            saddle_energy:  start_energy,
            barrier_energy: 0.0,
            start_energy,
            end_energy:     start_energy,
        };
        return Ok((
            vec![PathStep {
                structure:    start_dot,
                move_applied: None,
                energy:       start_energy,
                step_index:   0,
            }],
            stats,
        ));
    }

    // ── beam initialisation ───────────────────────────────────────────────────
    let mut beam = vec![PseudoIntermediate {
        pt:              start_pt.clone(),
        saddle_energy:   start_energy,
        remaining_moves: diff.move_list,
        path:            Vec::new(),
    }];

    // ── beam search ───────────────────────────────────────────────────────────
    if cache_stats {
        eprintln!(
            "# {:>4}  {:>7}  {:>7}  {:>6}  {:>7}  {:>5}  {:>6}  {:>6}",
            "step", "evals", "hits", "hit%", "unique", "surv", "errors", "pruned"
        );
    }
    let mut total_evals = 0usize;
    let mut total_hits  = 0usize;

    for step in 0..total_steps {
        // B2: collect lightweight expansions first; no path/remaining_moves clones yet.
        let mut expansions: Vec<Expansion> = Vec::new();

        let mut n_evals          = 0usize;
        let mut n_hits           = 0usize;
        let mut n_eval_errors    = 0usize;
        let mut n_ceiling_pruned = 0usize;

        for (parent_idx, parent) in beam.iter().enumerate() {
            for (move_idx, mv) in parent.remaining_moves.iter().enumerate() {
                let i = mv.i as usize;
                let j = mv.j as usize;

                // Validity: insertions require both positions unpaired;
                // deletions require the pair to exist.
                // Crossing pairs (pseudoknots) are allowed — no LoopTable check.
                if mv.is_insertion {
                    if parent.pt[i].is_some() || parent.pt[j].is_some() {
                        continue;
                    }
                } else if parent.pt[i] != Some(mv.j) {
                    continue;
                }

                // Apply move to a new PT.
                let mut new_pt = parent.pt.clone();
                if mv.is_insertion {
                    new_pt[i] = Some(mv.j);
                    new_pt[j] = Some(mv.i);
                } else {
                    new_pt[i] = None;
                    new_pt[j] = None;
                }

                // Memoized PK energy evaluation (B1: no string roundtrip).
                let (energy, hit) = match eval_energy(model, seq, &new_pt, &mut cache) {
                    Ok(pair) => pair,
                    Err(_) => { n_eval_errors += 1; continue; }
                };
                n_evals += 1;
                if hit { n_hits += 1; }

                // Energy ceiling filter.
                if max_energy.is_some_and(|cap| energy > cap) {
                    n_ceiling_pruned += 1;
                    continue;
                }

                let new_saddle = f64::max(parent.saddle_energy, energy);

                expansions.push(Expansion {
                    pt:             new_pt,
                    saddle_energy:  new_saddle,
                    current_energy: energy,
                    parent_idx,
                    move_idx,
                });
            }
        }

        total_evals += n_evals;
        total_hits  += n_hits;

        if expansions.is_empty() {
            let reason = if n_eval_errors > 0 && n_ceiling_pruned == 0 {
                format!(
                    "all {n_eval_errors} candidate move(s) produced intermediates the energy \
                     evaluator could not handle (complex pseudoknot topology with >2 bands or \
                     crossing strands). The target structure may require evaluator extensions."
                )
            } else if n_ceiling_pruned > 0 && n_eval_errors == 0 {
                format!(
                    "all {n_ceiling_pruned} candidate move(s) exceeded the energy ceiling \
                     ({:.2} kcal/mol). Try raising --max-energy or removing it entirely.",
                    max_energy.unwrap_or(f64::INFINITY)
                )
            } else {
                format!(
                    "no valid moves remain ({n_eval_errors} eval error(s), \
                     {n_ceiling_pruned} ceiling-pruned). \
                     Try removing --max-energy or checking the target topology."
                )
            };
            return Err(format!("Search stuck at step {step}: {reason}"));
        }

        // Sort: lowest saddle first; break ties by current energy.
        expansions.sort_by(|a, b| {
            a.saddle_energy
                .partial_cmp(&b.saddle_energy)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    a.current_energy
                        .partial_cmp(&b.current_energy)
                        .unwrap_or(Ordering::Equal)
                })
        });

        // Deduplicate by pair table (first occurrence = best path to that structure).
        // FxHashSet gives ~3× faster probes than SipHash for integer-heavy keys (B3).
        let mut seen: FxHashSet<PairTable> = FxHashSet::default();
        expansions.retain(|exp| seen.insert(exp.pt.clone()));

        // Prune to beam width.
        let n_unique    = expansions.len();
        expansions.truncate(beam_width);
        let n_survivors = expansions.len();

        if cache_stats {
            let hit_pct = if n_evals > 0 { 100.0 * n_hits as f64 / n_evals as f64 } else { 0.0 };
            eprintln!(
                "  {:>4}  {:>7}  {:>7}  {:>5.1}%  {:>7}  {:>5}  {:>6}  {:>6}",
                step + 1, n_evals, n_hits, hit_pct,
                n_unique, n_survivors, n_eval_errors, n_ceiling_pruned
            );
        }

        // B2: materialize only the survivors — at most `beam_width` clones.
        // B4: swap_remove(move_idx) is O(1) vs remove(move_idx) which is O(moves).
        beam = expansions.into_iter().map(|exp| {
            let parent = &beam[exp.parent_idx];
            let mv = parent.remaining_moves[exp.move_idx].clone();
            let mut new_remaining = parent.remaining_moves.clone();
            new_remaining.swap_remove(exp.move_idx);
            let mut new_path = parent.path.clone();
            new_path.push(mv);
            PseudoIntermediate {
                pt:              exp.pt,
                saddle_energy:   exp.saddle_energy,
                remaining_moves: new_remaining,
                path:            new_path,
            }
        }).collect();
    }

    // ── reconstruct full trajectory ───────────────────────────────────────────
    if cache_stats {
        let hit_pct = if total_evals > 0 { 100.0 * total_hits as f64 / total_evals as f64 } else { 0.0 };
        eprintln!(
            "# total: evals={total_evals} hits={total_hits} ({hit_pct:.1}%) cache_size={}",
            cache.len()
        );
    }

    let winner = beam.into_iter().next().ok_or("No beam survivors")?;
    reconstruct_path(model, seq, n, &start_pt, &start_dot, start_energy, &winner.path, &mut cache)
}

// ---------------------------------------------------------------------------
// Path reconstruction
// ---------------------------------------------------------------------------

fn reconstruct_path(
    model:        &ViennaRNA,
    seq:          &[ff_energy::Base],
    n:            usize,
    start_pt:     &PairTable,
    start_dot:    &str,
    start_energy: f64,
    moves:        &[Move],
    cache:        &mut FxHashMap<PairTable, f64>,
) -> Result<(Vec<PathStep>, PathStats), String> {
    let mut trajectory = Vec::with_capacity(moves.len() + 1);
    let mut pt = start_pt.clone();

    trajectory.push(PathStep {
        structure:    start_dot.to_owned(),
        move_applied: None,
        energy:       start_energy,
        step_index:   0,
    });

    let mut saddle = start_energy;

    for (step_idx, mv) in moves.iter().enumerate() {
        let i = mv.i as usize;
        let j = mv.j as usize;
        if mv.is_insertion {
            pt[i] = Some(mv.j);
            pt[j] = Some(mv.i);
        } else {
            pt[i] = None;
            pt[j] = None;
        }

        let (en, _) = eval_energy(model, seq, &pt, cache)?;
        if en > saddle { saddle = en; }

        let dot = pair_table_to_dot_bracket(&pt)
            .unwrap_or_else(|_| "?".repeat(n));

        trajectory.push(PathStep {
            structure:    dot,
            move_applied: Some(mv.clone()),
            energy:       en,
            step_index:   step_idx + 1,
        });
    }

    let stats = PathStats {
        saddle_energy:  saddle,
        barrier_energy: saddle - start_energy,
        start_energy,
        end_energy:     trajectory.last().unwrap().energy,
    };

    Ok((trajectory, stats))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::{ViennaRNA, parameters::{RNA_MT09, RNA_DP09}, NucleotideVec, parse_structure};

    fn model() -> ViennaRNA {
        ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09)
    }

    /// Empty target — trivial path, no steps.
    #[test]
    fn test_empty_target() {
        let m = model();
        let seq = "GCGAAACGC";
        let target = ".........";
        let (path, stats) = findpath_pseudo(&m, seq, None, target, 1, None, false).unwrap();
        assert_eq!(path.len(), 1);
        assert_eq!(stats.barrier_energy, 0.0);
    }

    /// Non-pseudoknotted target: result must match energy_of_pseudoknotted_structure.
    #[test]
    fn test_non_pk_target() {
        let m = model();
        let seq = "GGGAAACCC";
        let target = "(((...)))";
        let (path, _stats) = findpath_pseudo(&m, seq, None, target, 4, None, false).unwrap();

        assert_eq!(path.last().unwrap().structure, target);

        let seq_vec = NucleotideVec::try_from_rna(seq).unwrap();
        let loops = parse_structure(target).unwrap();
        let expected = m.energy_of_pseudoknotted_structure(&seq_vec, &loops).unwrap() as f64 / 100.0;
        assert!(
            (path.last().unwrap().energy - expected).abs() < 0.01,
            "energy mismatch: got {}, expected {expected}", path.last().unwrap().energy
        );
    }

    /// H-type pseudoknot target: path must reach target (compared by pair table).
    #[test]
    fn test_h_type_pk_target() {
        let m = model();
        let seq    = "GCGAUUUCUGACCGCUUUUUUGUCAG";
        let target = "[[[....(((((]]]......)))))";
        let (path, stats) = findpath_pseudo(&m, seq, None, target, 40, None, false).unwrap();

        // pair_table_to_dot_bracket may assign different bracket families than the input,
        // so compare pair tables rather than strings.
        let final_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(path.last().unwrap().structure.as_str()).unwrap()
        ).unwrap();
        let target_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(target).unwrap()
        ).unwrap();
        assert_eq!(final_pt, target_pt, "final pair table must match target");
        println!("H-type PK path ({} steps), saddle = {:.2} kcal/mol:",
            path.len() - 1, stats.saddle_energy);
        for step in &path {
            println!("  [{:2}] {} {:.2} kcal/mol",
                step.step_index, step.structure, step.energy);
        }
    }

    /// Non-empty start: fold from the second stem of the target PK (helix 2 preformed)
    /// into the full H-type pseudoknot by inserting the 3 pairs of helix 1.
    ///
    /// Because start is exactly helix 2 of the target and helix 1 pairs are
    /// simply inserted (no deletions), every intermediate is a valid H-type PK —
    /// the energy evaluation is guaranteed not to hit complex topologies.
    #[test]
    fn test_non_empty_start() {
        let m = model();
        // Sequence: GCGAUUUCUGACCGCUUUUUUGUCAG  (26 nt)
        // Start:    .......(((((.........)))))  — helix 2 only, pairs (7,25)...(11,21)  [26 nt]
        // Target:   [[[....(((((]]]......)))))  — full H-type PK
        // Diff: 3 insertions (0,12),(1,13),(2,14); 0 deletions
        let seq    = "GCGAUUUCUGACCGCUUUUUUGUCAG";
        let start  = ".......(((((.........)))))";
        let target = "[[[....(((((]]]......)))))";
        let (path, stats) = findpath_pseudo(&m, seq, Some(start), target, 10, None, false).unwrap();

        // First step must be the start structure (by pair table).
        let first_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(path.first().unwrap().structure.as_str()).unwrap()
        ).unwrap();
        let start_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(start).unwrap()
        ).unwrap();
        assert_eq!(first_pt, start_pt, "first step must be the start structure");

        // Last step must be the target structure.
        let final_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(path.last().unwrap().structure.as_str()).unwrap()
        ).unwrap();
        let target_pt = extended_dot_bracket_to_pair_table(
            &ExtendedDotBracketVec::try_from(target).unwrap()
        ).unwrap();
        assert_eq!(final_pt, target_pt, "final pair table must match target");

        println!("Non-empty start PK path ({} steps), saddle = {:.2} kcal/mol:",
            path.len() - 1, stats.saddle_energy);
        for step in &path {
            println!("  [{:2}] {} {:.2} kcal/mol",
                step.step_index, step.structure, step.energy);
        }
    }
}
