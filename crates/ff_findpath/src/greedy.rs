use crate::utils::{Intermediate, PathStats, PathStep, compare_structures, apply_move};
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel};
use ff_structure::{PairTable, DotBracketVec};

/// Calculates a folding path and an energy barrier between two RNA secondary structures using the greedy heuristics (Morgan-Higgs algorithm).
/// # Input 
/// It takes a 
/// -  `model`: The ViennaRNA energy model to evaluate energies., 
/// -  `sequence`: the RNA sequence as a string,
/// -  `s1`: the starting structure in dot-bracket notation, 
/// -  `s2`: the target structure in dot-bracket notation,

/// 
/// # Output
/// It returns a `Result` which is either:
/// - `Ok((Vec<PathStep>, PathStats))` containing the folding path (as a vector of `PathStep`) and some statistics about the path (as `PathStats`), or 
/// - `Err(String)` containing an error message if the input parameters are invalid or if the search fails (e.g. due to energy constraints or topological issues).
/// 
/// # Internal calls:
/// - `compare_structures()`: This function compares the starting and target structures to generate a list of moves that can transform one structure into the other. This is done once at the beginning to prepare the move lists for both directions.
/// - `apply_move()`: This function generates neighboring structures by applying the allowed moves to the current structure, while also calculating their energies
/// 
/// NOTE: the `max_energy` parameter is not used in greedy search, as it is a greedy algorithm that explores all valid neighbors and picks the best one.
pub fn greedy_find_path(
    model: &ViennaRNA,
    sequence: &str,
    s1: &str,
    s2: &str,
) -> Result<(Vec<PathStep>, PathStats), String> {
    
    // 1. Setup & Validation
    let seq_vec = NucleotideVec::try_from_rna(sequence).expect("Failed to parse RNA sequence");
    let pt_start = PairTable::try_from(s1).map_err(|_| "Invalid start structure s1")?;
    let pt_target = PairTable::try_from(s2).map_err(|_| "Invalid target structure s2")?;

    // 2. Prepare Moves
    let comparison = compare_structures(&pt_start, &pt_target);
    let move_list = comparison.move_list;
    let total_steps = move_list.len();

    // 3. Initialize Energies
    let start_energy = model.energy_of_structure(&seq_vec, &pt_start).expect("failed to calculate starting energy") as f64 / 100.0;
    
    // 4. Initialize Trajectory
    let mut trajectory = Vec::with_capacity(total_steps + 1);
    
    // Record step 0 (Start)
    trajectory.push(PathStep {
        structure: DotBracketVec::try_from(&pt_start).unwrap().to_string(), 
        move_applied: None,
        energy: start_energy,
        step_index: 0,
    });

    // 5. Initialize Search State
    let mut current = Intermediate {
        pt: pt_start,
        saddle_energy: start_energy,
        current_energy: start_energy,
        remaining_moves: move_list,
        path: Vec::new(),
    };

    // 6. Greedy Loop
    for step in 1..=total_steps {
        
        // A. Generate Neighbors
        let candidates = apply_move(model, &seq_vec, &current, None);

        if candidates.is_empty() {
            return Err(format!("Greedy search stuck at step {}: No valid topology neighbors.", step));
        }

        // B. Select Best Neighbor (Min saddle, then Min current)
        let best_candidate = candidates
            .into_iter()
            .min_by(|a, b| {
                a.saddle_energy.partial_cmp(&b.saddle_energy)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| {
                        a.current_energy.partial_cmp(&b.current_energy)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            })
            .unwrap();

        // C. Record Step
        let last_move = best_candidate.path.last().cloned();
        
        trajectory.push(PathStep {
            structure: DotBracketVec::try_from(&best_candidate.pt).unwrap().to_string(),
            move_applied: last_move,
            energy: best_candidate.current_energy,
            step_index: step,
        });

        // D. Advance State
        current = best_candidate;
    }

    // 7. Finalize Stats
    let stats = PathStats {
        saddle_energy: current.saddle_energy,
        barrier_energy: current.saddle_energy - start_energy,
        start_energy: start_energy,
        end_energy: current.current_energy,
    };

    Ok((trajectory, stats))
}




// #####################################
//                  TESTS
// #####################################




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_greedy_find_path_integration_01() {
        // Simple hairpin closing example
        // Assume ViennaRNA default model
        
        let seq = "GGGGUUUUCCCC";
        let s1 = "............"; //(Unpaired)
        let s2 = "((((....))))"; //(Full Hairpin)
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();
        let result = greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((trajectory, stats)) => {
                println!("Path found with {} steps", trajectory.len());
                println!("Saddle Energy: {:.2} kcal/mol", stats.saddle_energy);
                println!("Barrier Height: {:.2} kcal/mol", stats.barrier_energy);
                
                // Assertions
                assert_eq!(trajectory.len(), 5, "Should have 5 PathStep entries in the trajectory (including start and end)");
                // Check if the final path actually formed the structure
                // We know it starts empty, so all moves should be insertions
                
                for m in trajectory{
                    println!("{:?}", m)
                }
                }
            
            Err(e) => panic!("Greedy path search failed: {}", e),
        }
    }

    #[test]
    fn test_greedy_find_path_integration_02() {

       // Assume ViennaRNA default model
        
        let seq = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let s1 = ".(((..............)))....";
        let s2 = "..((((.........))))......";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();
        let result = greedy_find_path(&model, seq, s1, s2);

        match result {
            Ok((trajectory, stats)) => {
                println!("Path found with {} steps", trajectory.len());
                println!("Saddle Energy: {:.2} kcal/mol", stats.saddle_energy);
                println!("Barrier Height: {:.2} kcal/mol", stats.barrier_energy);
                
                // Assertions
                assert_eq!(trajectory.len(), 8, "Should have 8 PathStep entries in the trajectrory (including start and end)");
                // Check if the final path actually formed the structure
                for m in trajectory{
                    println!("{:?}", m)
                }
                
            },
            Err(e) => panic!("Greedy path search failed: {}", e),
        }
    }

    #[test]
    fn test_greedy_find_path_printing_debugging(){
        // Does not test anything, just prints the path and stats for visual inspection and debugging.
        let seq = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let s1 = ".(((..............)))....";
        let s2 = "..((((.........))))......";
        
        // Initialize Model (This might fail if params aren't found in your env)
        // Ensure you have valid params or a mock model available.
        let model = ViennaRNA::default();

        let (steps, stats) = greedy_find_path(&model, seq, s1, s2).unwrap();
        println!("Test steps:");
        for step in steps {
            println!("{} \t {} \t {}", step.structure,step.move_applied.unwrap_or_default(), step.energy );
        } 
        println!("Stats: {}", stats);
    }
}


