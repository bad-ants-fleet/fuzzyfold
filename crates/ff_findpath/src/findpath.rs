use crate::utils::{Move, Intermediate, PathStats, PathStep, compare_structures, apply_move};
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel};
use ff_structure::{PairTable, DotBracketVec};
use std::cmp::Ordering;
use std::collections::HashSet;


//NOTE: the apply_move function in greedy.rs is implemented to have a maximum energy barrier parameter, which is used in findpath but not in the greedy algorithm.

/// Calculates a folding path and an energy barrier between two RNA secondary structures using the findpath algorithm.
/// # Input 
/// It takes a 
/// -  `model`: The ViennaRNA energy model to evaluate energies., 
/// -  `sequence`: the RNA sequence as a string,
/// -  `s1`: the starting structure in dot-bracket notation, 
/// -  `s2`: the target structure in dot-bracket notation,
/// -  `target_m`: the desired search width parameter `m` for the findpath algorithm (default: 1 = greedy search), 
/// -  `max_energy`: an optional parameter to set a maximum energy threshold for the search. 
/// 
/// 
/// # Output
/// It returns a `Result` which is either:
/// - `Ok((Vec<PathStep>, PathStats))` containing the folding path (as a vector of `PathStep`) and some statistics about the path (as `PathStats`), or 
/// - `Err(String)` containing an error message if the input parameters are invalid or if the search fails (e.g. due to energy constraints or topological issues).
/// 
/// 
/// # Example
/// ```rust
/// use ff_findpath::findpath::findpath;
/// use ff_energy::ViennaRNA;
/// 
/// let model = ViennaRNA::default();
/// let sequence = "AGCCAUGAGUGUAUAGUGGGCCUAU";
/// let struct1 = ".(((..............)))....";
/// let struct2 = "..((((.........))))......";
/// let result = findpath(&model, sequence, struct1, struct2, 10, None);
/// 
/// let (path, stats) = result.unwrap();
/// println!("Folding Path found with findpath algorithm:");
/// println!("-----------------");
/// println!("{} \t applied move \t energy", sequence);
/// for step in path {
///     println!("{} \t {} \t {} kcal/mol", step.structure, step.move_applied.unwrap_or_default(), step.energy );
/// }   
/// println!("-----------------");
/// println!("Statistics:");
/// println!("{}", stats);
/// ```
/// 
/// # Algorithmic details
/// The search-with parameter `m` is used to increase the search space by allowing the `m` best neighbors (in terms of lowest energy barrier) to be considered at each step, instead of just the single best neighbor (i.e `m=1` is equivalent to the greedy algorithm). 
/// The algorithm iteratively expands the search space by increasing `m` until it reaches the `target_m` specified by the user. 
/// The search space is explored using a **breath-first-search** approach, where the search-with parameter is doubled at each iteration to allow for wider search.
///
/// It starts with with `m=1` (greedy search) and doubles `m` in each iteration, allowing for a wider search and potentially upper energy barriers to be overcome.
/// This can help to find better paths that might be missed by a purely greedy approach, at the cost of increased computational time.
/// Also the search is performed from both directions (s1 to s2 and s2 to s1), alternating after each applied move (bidirectional beam search).
/// The energy barrier is tracked during the search, and if a candidate move results in a structure with energy above the specified `max_energy` threshold, that move is discarded and not added to the beam for the next iteration.
/// Initally the `max_energy` threshold is not set, but after the first successful pass, it is updated to the saddle energy of the best path found so far. 
///
///
/// # Internal function calls 
/// - `compare_structures()`: This function compares the starting and target structures to generate a list of moves that can transform one structure into the other. This is done once at the beginning to prepare the move lists for both directions.
/// - `run_beam_pass()`: This function runs a single, direction-agnostic beam search pass with a specfiic `m` value. 
/// - `invert_path_trajectory()`: This function takes a backward trajectory (s2 -> s1), extracts the moves, reverses and inverts them, and reconstructs the forward trajectory (s1 -> s2).
pub fn findpath(
    model: &ViennaRNA,
    sequence: &str,
    s1: &str,
    s2: &str,
    target_m: usize,
    mut max_energy: Option<f64>, 
) -> Result<(Vec<PathStep>, PathStats), String> {
    
    // 1. Setup & Pre-calculation (Done ONCE)
    let seq_vec = NucleotideVec::try_from_rna(sequence)
        .expect("Failed to parse RNA sequence");
    
    let pt_start = PairTable::try_from(s1).map_err(|_| "Invalid start structure s1")?;
    let pt_target = PairTable::try_from(s2).map_err(|_| "Invalid target structure s2")?;

    // Calculate Moves ONCE (Forward: s1 -> s2)
    let fwd_comp = compare_structures(&pt_start, &pt_target);
    let moves_fwd = fwd_comp.move_list;

    // Derive Backward Moves (s2 -> s1) by inverting the insertion logic
    let moves_bwd: Vec<Move> = moves_fwd.iter().map(|m| Move {
        i: m.i,
        j: m.j,
        is_insertion: !m.is_insertion,
    }).collect();

    // 2. State Variables
    let mut current_m: usize = 1;
    let mut forward_dir = true; // True = s1->s2, False = s2->s1
    // We add a bool to explicitly store whether this specific result was forward or backward (final reconstruction needs to know this)
    let mut last_result: Option<(Vec<PathStep>, PathStats, bool)> = None;

    // 3. Iterative Loop
    // start with m = 1, and keep increasing until we reach target_m.
    // also, we keep alternating between forward and backward search.
    loop {
        // A. Select Inputs based on Direction
        let (current_pt, current_moves) = if forward_dir {
            (&pt_start, &moves_fwd)
        } else {
            (&pt_target, &moves_bwd)
        };

        // B. Run the Direction-Agnostic Pass
        let pass_result = run_beam_pass(
            model, 
            &seq_vec, 
            current_pt, 
            current_moves, 
            current_m, 
            max_energy
        );
        // match on the results of the last pass. 
        // If it was successful, we update max_energy and save the path + stats + direction. 
        // If it failed, we do nothing and let the loop continue to the next iteration (with increased m and flipped direction). The last successful result remains safely stored.
        match pass_result {
            Ok((path, stats)) => {
                max_energy = Some(stats.saddle_energy);
                // SAVE the direction (forward_dir) along with the path!
                last_result = Some((path, stats, forward_dir));
            },
            Err(_e) => {
                // TODO ask Stefan what to do in case of a failed pass (e.g. due to energy constraint or topology.
                // Do nothing! Let the loop continue to the next width (m) 
                // and direction. The last successful result remains safely stored.
            }
        }
            
            


        // C. Check Exit Condition
       if current_m >= target_m { 
            break;
        }

        // D. Prepare Next Iteration
        // Double m, set m to target_m if doubling would exceed it.
        let next_m = current_m * 2;
        current_m = if next_m > target_m { target_m } else { next_m };

        // Flip Direction
        forward_dir = !forward_dir;

    }

    // 4. Finalize Result
    // Unpack the saved `was_fwd` flag here
    if let Some((path, stats, was_fwd)) = last_result {
        if was_fwd {
            // Path was found S1 -> S2
            Ok((path, stats))
        } else {
            // Path was found S2 -> S1, we must invert it!
            invert_path_trajectory(model, &seq_vec, &pt_start, &path)
        }
    } else {
        Err("Search loop finished without producing a result. All passes failed.".to_string())
    }
}



// ------------ Worker Function ------------
/// Runs a single, direction-agnostic beam search pass.
/// 
/// # Input:
/// - `model`: The ViennaRNA energy model to evaluate energies.
/// - `seq_vec`: The sequence as a vector of nucleotides.
/// - `start_pt`: The starting structure as a PairTable.
/// - `initial_moves`: The list of moves to apply (direction-agnostic).
/// - `m`: The beam width for this pass.
/// - `max_energy`: Optional maximum energy threshold to filter candidates.
/// 
/// # Output:
/// - `Ok((Vec<PathStep>, PathStats))` if a path is found within the constraints.
/// - `Err(String)` if the search fails (e.g. due to energy constraints or topological issues).
/// 
/// # Algorithmic details:
/// The function initializes a beam with the starting structure and iteratively expands it by applying the allowed moves. 
/// At each step, it generates candidate neighbors, filters them based on the `max_energy`
/// threshold, sorts them by saddle energy (and current energy as a tiebreaker), deduplicates them, and prunes to keep only the top `m` candidates for the next iteration.
/// 
/// # Internal function calls:
/// - `apply_move()`: This function generates neighboring structures by applying the allowed moves to the current structure, while also calculating their energies and filtering based on the `max_energy` threshold.
/// - `reconstruct_path_from_moves()`: If a valid path is found after the final iteration, this function reconstructs the full trajectory of structures and energies from the sequence of moves that led to the solution.
pub fn run_beam_pass(
    model: &ViennaRNA,
    seq_vec: &NucleotideVec,
    start_pt: &PairTable,   
    initial_moves: &[Move], 
    m: usize,
    max_energy: Option<f64>,
) -> Result<(Vec<PathStep>, PathStats), String> {

    let total_steps = initial_moves.len();
    let start_energy = model.energy_of_structure(seq_vec, start_pt).expect("Failed to calculate starting energy") as f64 / 100.0;

    // Check if the start structure already violates max_energy
    if let Some(max_e) = max_energy {
        if start_energy > max_e {
            return Err(format!("Start energy {:.2} already exceeds max_energy {:.2}", start_energy, max_e));
        }
    }

    let mut beam = vec![Intermediate {
        pt: start_pt.clone(),
        saddle_energy: start_energy,
        current_energy: start_energy,
        remaining_moves: initial_moves.to_vec(), 
        path: Vec::new(),
    }];

    for _step in 1..=total_steps {
        let mut next_candidates = Vec::new();

        // 1. Expand all survivors in the beam
        for parent in beam {
            let mut neighbors = apply_move(model, seq_vec, &parent, max_energy);
            next_candidates.append(&mut neighbors);
        }

        if next_candidates.is_empty() {
             return Err("Search stuck (dead end due to topology or energy constraint)".to_string());
        }

        // 2. Sort: Min Saddle -> Min Current
        // This ensures the BEST paths are at the front of the list.
        next_candidates.sort_by(|a, b| {
            a.saddle_energy.partial_cmp(&b.saddle_energy).unwrap_or(Ordering::Equal)
                .then_with(|| a.current_energy.partial_cmp(&b.current_energy).unwrap_or(Ordering::Equal))
        });

        // 3. Deduplicate
        // We track the PairTables we've seen in this specific generation.
        // Because the list is sorted, the first time we `insert` a structure, 
        // it is the lowest-energy path to that structure.
        let mut seen_structures = HashSet::new();
        
        next_candidates.retain(|candidate| {
            // Note: This assumes PairTable implements `Eq` and `Hash`.
            // If it doesn't, you can hash a string representation:
            // let hash_key = DotBracketVec::try_from(&candidate.pt).unwrap().to_string();
            // seen_structures.insert(hash_key)
            
            seen_structures.insert(candidate.pt.clone()) 
        });

        // 4. Prune to keep only top 'm' unique candidates
        if next_candidates.len() > m {
            next_candidates.truncate(m);
        }

        beam = next_candidates;
    }

    // Select the absolute best winner
    let winner = beam.first().ok_or("No survivors in beam")?;
    
    // Reconstruct full trajectory history
    reconstruct_path_from_moves(model, seq_vec, start_pt, &winner.path)
}





// ----------------- Helpers ----------------

/// Helper function:
/// Reconstructs the full trajectory from a sequence of moves.
/// 
/// # Input:
/// - `model`: The ViennaRNA energy model to evaluate energies.
/// - `seq_vec`: The sequence as a vector of nucleotides.
/// - `start_pt`: The starting structure as a PairTable.
/// - `moves`: The sequence of moves that led to the solution.
/// 
/// # Output:
/// - `Ok((Vec<PathStep>, PathStats))` containing the full trajectory and statistics.
/// - `Err(String)` if reconstruction fails (e.g. due to invalid moves).
///
pub fn reconstruct_path_from_moves(
    model: &ViennaRNA,
    seq_vec: &NucleotideVec,
    start_pt: &PairTable,
    moves: &[Move],
) -> Result<(Vec<PathStep>, PathStats), String> {
    
    let mut trajectory = Vec::with_capacity(moves.len() + 1);
    let mut pt = start_pt.clone();
    let start_en = model.energy_of_structure(seq_vec, &pt).expect("Failed to calculate starting energy") as f64 / 100.0;

    trajectory.push(PathStep {
        structure: DotBracketVec::try_from(&pt).unwrap().to_string(), 
        move_applied: None,
        energy: start_en,
        step_index: 0,
    });

    let mut saddle = start_en;

    for (i, mv) in moves.iter().enumerate() {
        // Apply move logic
        if mv.is_insertion {
            pt[mv.i] = Some(mv.j as u16);
            pt[mv.j] = Some(mv.i as u16);
        } else {
            pt[mv.i] = None;
            pt[mv.j] = None;
        }

        // Calculate and track energy
        let en = model.energy_of_structure(seq_vec, &pt).expect("Failed to calculate energy") as f64 / 100.0;
        if en > saddle { saddle = en; }

        trajectory.push(PathStep {
            structure: DotBracketVec::try_from(&pt).unwrap().to_string(), 
            move_applied: Some(mv.clone()),
            energy: en,
            step_index: i + 1,
        });
    }

    let stats = PathStats {
        saddle_energy: saddle,
        barrier_energy: saddle - start_en,
        start_energy: start_en,
        end_energy: trajectory.last().unwrap().energy,
    };

    Ok((trajectory, stats))
}

/// Helper function:
/// Takes a backward trajectory (s2 -> s1) and inverts it to reconstruct the forward trajectory (s1 -> s2).
/// 
/// # Input:
/// - `model`: The ViennaRNA energy model to evaluate energies.
/// - `seq_vec`: The sequence as a vector of nucleotides.
/// - `true_start_pt`: The true starting structure (s1) as a PairTable
/// - `wrong_dir_path`: The trajectory found in the wrong direction (s2 -> s1) as a slice of PathStep.
/// 
/// # Output:
/// - `Ok((Vec<PathStep>, PathStats))` containing the reconstructed forward trajectory and statistics.
/// - `Err(String)` if the inversion or reconstruction fails (e.g. due to invalid moves).
/// 
/// # Internal function calls:
/// - `reconstruct_path_from_moves()`: reconstruct the forward trajectory
pub fn invert_path_trajectory(
    model: &ViennaRNA,
    seq_vec: &NucleotideVec,
    true_start_pt: &PairTable,
    wrong_dir_path: &[PathStep],
) -> Result<(Vec<PathStep>, PathStats), String> {
    
    // Extract moves (skipping the step 0 start state)
    let moves: Vec<Move> = wrong_dir_path.iter()
        .skip(1) 
        .filter_map(|step| step.move_applied.clone())
        .collect();

    // Reverse the order and invert the insertion/deletion logic
    let inverted_moves: Vec<Move> = moves.iter().rev().map(|m| Move {
        i: m.i, 
        j: m.j, 
        is_insertion: !m.is_insertion
    }).collect();

    // Reconstruct the proper trajectory starting from S1
    reconstruct_path_from_moves(model, seq_vec, true_start_pt, &inverted_moves)
}








// #####################################
//                  TESTS
// #####################################







#[cfg(test)]
mod tests {
    //use super::findpath;
    use super::findpath;
    use ff_energy::ViennaRNA;

    #[test]
    fn test_findpath_basic() {
        let model = ViennaRNA::default();
        let sequence = "AGCCAUGAGUGUAUAGUGGGCCUAU"; // short.txt from test_data
        let s1 = ".(((..............)))...."; 
        let s2 = "..((((.........))))......"; 

        let result = findpath(&model, sequence, s1, s2, 1, None);
        assert!(result.is_ok(), "Expected a valid path, got error: {:?}", result.err());
        
        let (path, stats) = result.unwrap();
        assert_eq!(path.len(), 8); // 
        assert_eq!(stats.start_energy, -2.2);
        assert_eq!(stats.saddle_energy, 4.3);
        println!("Test steps:");
        
        for step in path {
            println!("{} \t {} \t {}", step.structure,step.move_applied.unwrap_or_default(), step.energy );
        } 
        println!("Stats: {:?}", stats);
    }
}