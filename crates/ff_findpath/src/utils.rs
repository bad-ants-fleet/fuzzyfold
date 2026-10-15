use std::fmt;
use ff_energy::{NucleotideVec, ViennaRNA, EnergyModel}; 
use ff_structure::{LoopInfo, LoopTable, NAIDX, PairTable};

/// Corresponds to an elementary move in the transformation from two RNA structures.
///
/// `Move` represents the insertion or deletion of a base pair at specific indices (i, j) between two structures. 
/// 
/// # Fields:
/// - `i`: The index of the first base in the pair.
/// - `j`: The index of the second base in the pair.
/// - `is_insertion`: A boolean flag indicating whether the move is an insertion (true) or deletion (false).
/// 
/// 
/// For example, if we have a move with i = 2, j = 5, is_insertion = true, it means that we are inserting the base pair (2, 5) into the structure.
/// The corresponding PairTable representation would have `pt[2] = Some(5)` and `pt[5] = Some(2)`.
/// 
/// # Note:
/// We define the base pair distance between two structures as the number of moves required to transform one structure into the other.
/// This means one move corresponds to the change of two indices in the PairTable representation.
/// 
/// # Example usage:
/// ```rust
/// use ff_findpath::utils::Move;
/// 
/// let mv = Move { i: 2, j: 5, is_insertion: true };
/// ```
/// 
// Derive Default, makes default Move { i: 0, j: 0, is_insertion: false } (del (0,0) at start of every path)
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct Move {
    pub i: NAIDX,
    pub j: NAIDX,
    pub is_insertion: bool, // true = insert, false = delete
}

// Display implementations for better readability
impl fmt::Display for Move {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = if self.is_insertion { "ins" } else { "del" };
        write!(f, "{}({}, {})", kind, self.i, self.j)
    }
}



/// StructureDifference captures the difference between two secondary structures.
/// 
/// # Fields:
/// - `move_list`: A list of moves (insertions or deletions of base pairs) required to transform one structure into the other.
/// - `bp_distance``: The base pair distance between the two structures, defined as the number of base pairs that differ between them.
/// 
/// # Note:
/// Each move in the move_list corresponds to a change of two indices in the PairTable representation,
/// hence a bp_distance of one equals one move.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureDifference {
    pub move_list: Vec<Move>,
    pub bp_distance: u32,
}

impl fmt::Display for StructureDifference {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "bp_distance: {}", self.bp_distance)?;
        // don't print the number of moves because they are equal to bp_distance
        writeln!(f, "moves ({}):", self.move_list.len())?;
        for mv in &self.move_list {
            writeln!(f, "  - {}", mv)?;
        }
        Ok(())
    }
}


/// Represents the current structure and state of the search during a folding path calculation.
/// 
/// This struct is used to keep track of the current structure (as a PairTable), the saddle energy encountered so far, the current energy, the list of remaining moves that can be applied, and the path of moves taken to reach this state.
/// # Fields:
/// - `pt`: The current structure represented as a PairTable.
/// - `saddle_energy`: The highest energy encountered along the path to reach this structure (used for calculating the energy barrier).
/// - `current_energy`: The energy of the current structure.
/// - `remaining_moves`: A list of moves that can still be applied to transform the current structure towards the target structure. 
/// - `path`: The sequence of moves that have been applied to reach the current structure from the starting structure. This is used for reconstructing the folding path once the target structure is reached.
/// 
/// 
#[derive(Clone, Debug)]
pub struct Intermediate {
    pub pt: PairTable,
    pub saddle_energy: f64,
    pub current_energy: f64,
    pub remaining_moves: Vec<Move>, // Moves available to be taken
    pub path: Vec<Move>,            // Sequence of moves taken so far
}



/// Represents a single step in the folding trajectory for analysis.
/// 
/// Each `PathStep` captures the structure at that step, the move that was applied to get there (if any), and the energy of that structure.
#[derive(Debug, Clone)]
pub struct PathStep {
    pub structure: String,      // Dot-bracket representation
    pub move_applied: Option<Move>, // The move that led to this state (None for start)
    pub energy: f64,            // Energy of this structure
    pub step_index: usize,
}

impl fmt::Display for PathStep {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let move_str = match &self.move_applied {
            Some(mv) => format!("{}", mv),
            None => "Start".to_string(),
        };
        write!(
            f,
            "Step {}: Structure: {}, Move: {}, Energy: {:.2}",
            self.step_index, self.structure, move_str, self.energy
        )
    }
}

/// Summary statistics for the entire path.
/// 
/// Includes the saddle energy (highest energy point), barrier energy (saddle - start), start energy, and end energy.
#[derive(Debug, Clone)]
pub struct PathStats {
    pub saddle_energy: f64,     // Highest energy point (max_en)
    pub barrier_energy: f64,    // max_en - start_en
    pub start_energy: f64,
    pub end_energy: f64,
}

impl fmt::Display for PathStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Saddle energy: {:.2} kcal/mol, Barrier energy: {:.2} kcal/mol, Start energy: {:.2} kcal/mol, End energy: {:.2} kcal/mol",
            self.saddle_energy,
            self.barrier_energy,
            self.start_energy,
            self.end_energy
        )
    }
}




/// A lightweight struct to hold a potential neighbor of a secondary structure.
/// 
/// This struct is used during the generation of neighbors in the search algorithms. It captures the new PairTable resulting from applying a move, the move that was applied, and the original index of that move in the list of available moves (which is useful for tracking and pruning).
pub struct NeighborCandidate {
    pub pt: PairTable,
    pub applied_move: Move,
    pub original_move_index: usize, // Needed to remove from remaining_moves later
}




/// Calculates the difference between two structures and prepares the list of moves required to transform structure 1 into structure 2
/// 
/// # Input
/// - `pt1`: The first structure represented as a PairTable. This is the structure from which we want to transform.
/// - `pt2`: The second structure represented as a PairTable. This is the structure to which we want to transform.
/// 
/// # Output
/// It returns a `StructureDifference` struct containing the list of moves required to transform pt1 into pt2, as well as the base pair distance between the two structures.
pub fn compare_structures(pt1: &PairTable, pt2: &PairTable) -> StructureDifference {
    let mut diff = StructureDifference {
        move_list: vec![],
        bp_distance: 0,
    };

    let length = pt1.len(); // Both structures should have the same length, we can check this in the caller function if needed

    // Iterate through all bases
    for i in 0..length {
        // We only care if the pairing status at i is different
        if pt1[i] != pt2[i] {
            
            // 1. Check for Deletion from pt1
            // logic: if pt1 has a pair (i, j) starting at i
            if let Some(j) = pt1[i] {
                // Ensure we only process the pair once (at the opening bracket, where i < j)
                if i < j as usize {
                    let deletion_move = Move {
                        i: i as NAIDX,
                        j: j as NAIDX,
                        is_insertion: false,
                    };

                    diff.move_list.push(deletion_move);
                    diff.bp_distance += 1;
                }
            }

            // 2. Check for Insertion into pt2
            // logic: if pt2 has a pair (i, j) starting at i
            if let Some(j) = pt2[i] {
                // Ensure we only process the pair once (at the opening bracket)
                if i < j as usize {
                    let insertion_move = Move {
                        i: i as NAIDX,
                        j: j as NAIDX,
                        is_insertion: true,
                    };

                    
                    diff.move_list.push(insertion_move);
                    diff.bp_distance += 1;
                }
            }
        }
    }

    // Python does `move_list.sort(key=lambda x: x[0])`
    // We sort by the index `i` to maintain consistent order
    // TODO: is sorting useful, because if a insertion of a base pair (i, j) is before a deletion that affects the same indices,
    // then the order matters? ==> maybe other function solves this.
    
    //diff.move_list.sort_by_key(|m| m.i);
    // 2. Sort moves directly
    // Primary Key: Deletion (false) before Insertion (true)
    // Secondary Key: Index i (ascending)
    // Tertiary Key: Index j (ascending)
    diff.move_list.sort_by(|a, b| {
        a.is_insertion.cmp(&b.is_insertion)
            .then(a.i.cmp(&b.i))
            .then(a.j.cmp(&b.j))
    });
    diff
}

/// Wrapper function to get the list of moves required to transform one structure into the other
/// 
/// Simply calls the `compare_structures` function and extracts the move_list from the returned `StructureDifference`.
/// 
/// # Input
/// - `pt1`: The first structure represented as a PairTable. This is the structure from which we want to transform.
/// - `pt2`: The second structure represented as a PairTable. This is the structure to which we want to transform.
/// 
#[allow(dead_code)]
fn prepare_moves(pt_1: &PairTable, pt_2: &PairTable) -> Vec<Move> {
    // compare_structures returns StructureDifference { move_list, bp_distance }
    // We only want the move_list 
    let diff = compare_structures(pt_1, pt_2);
    
    diff.move_list
}








/// Generates all valid neighbors of the given input structure.
/// 
/// Checks if the move is valid (topology) and applies it to generate a new PairTable, as well as the move that was valid.
/// 
/// NOTE: This function does NOT calculate energy or clone paths. It only checks if the move is valid (e.g. no pseudoknots) and applies it to generate a new PairTable. 
/// The energy calculation and path cloning are done in the apply_move function, which calls this function to get valid neighbors before calculating their energies and updating the search state.
/// 
/// # Input
/// - `current_pt`: The current structure represented as a PairTable. This is the structure from which we want to generate neighbors by applying the remaining moves.
/// - `available_moves`: A list of moves that can potentially be applied to the current structure to generate neighbors. 
/// 
/// # Output
/// It returns a vector of `NeighborCandidate`, where each candidate includes the new PairTable resulting from applying a valid move, the move that was applied, and the original index of that move in the `available_moves` list.
pub fn generate_valid_neighbors(
    current_pt: &PairTable,
    available_moves: &[Move],
) -> Vec<NeighborCandidate> {
    let mut candidates = Vec::new();

    // LoopTable is calculated once per expansion (O(N))
    let loop_table = match LoopTable::try_from(current_pt) {
        Ok(lt) => lt,
        Err(_) => return Vec::new(), // Invalid current structure, no neighbors
    };

    for (idx, candidate_move) in available_moves.iter().enumerate() {
        let i = candidate_move.i;
        let j = candidate_move.j;

        // --- VALIDITY CHECKS (Topology) ---
        if candidate_move.is_insertion {
            // Check for pseudoknots: i and j must be in the same loop context
            match (&loop_table[i as usize], &loop_table[j as usize]) {
                (LoopInfo::Unpaired { l: l1 }, LoopInfo::Unpaired { l: l2 }) => {
                    if l1 != l2 { continue; } // Different loops = crossing = invalid
                },
                _ => continue, // Already paired
            }
        } else {
            // Deletion: The pair must actually exist
            if current_pt[i as usize] != Some(j) { continue; }
        }

        // --- APPLY MOVE ---
        let mut new_pt = current_pt.clone();
        if candidate_move.is_insertion {
            new_pt[i as usize] = Some(j);
            new_pt[j as usize] = Some(i);
        } else {
            new_pt[i as usize] = None;
            new_pt[j as usize] = None;
        }

        candidates.push(NeighborCandidate {
            pt: new_pt,
            applied_move: candidate_move.clone(),
            original_move_index: idx,
        });
    }

    candidates
}

/// Performs one expansion step from the given Intermediate state.
/// 
/// It generates all valid neighbors by calling the `generate_valid_neighbors` function, calculates their energies,
///
/// # Input
/// - `model`: The ViennaRNA energy model to evaluate energies.
/// - `seq_vec`: The RNA sequence as a vector of nucleotides (NucleotideVec).
/// - `intermediate`: The current state of the search, including the current structure (as a PairTable), the saddle energy so far, the current energy, the list of remaining moves, and the path of moves taken to reach this state.
/// - `max_energy`: An optional parameter to set a maximum energy threshold for the search. If provided, any neighbor with energy above this threshold will be filtered out and not returned as a valid neighbor. 
/// This can be used to prune the search space and focus on more promising trajectories, but it is not used in the greedy algorithm (as it explores all valid neighbors and picks the best one).
/// 
/// # Output
/// It returns a vector of `Intermediate` states, each representing a valid neighbor of the input `intermediate`.
/// 
/// # Internal calls:
/// - generate_valid_neighbors to get valid topologies
pub fn apply_move(
    model: &ViennaRNA,
    seq_vec: &NucleotideVec,
    intermediate: &Intermediate,
    max_energy: Option<f64>,
) -> Vec<Intermediate> {
    
    // 1. Get valid topologies
    let candidates = generate_valid_neighbors(&intermediate.pt, &intermediate.remaining_moves);
    let mut results = Vec::new();

    for cand in candidates {
        // 2. Calculate Energy (Expensive)
        let en = model.energy_of_structure(seq_vec, &cand.pt).expect("failed to calculate Energy") as f64 / 100.0;
        

        // NOT USED IN GREEDY as we always explore all moves and pick the best
        // 3. Filter (Optimization: Drop before cloning paths if eneregy too high)
        if let Some(max) = max_energy {
            if en >= max { continue; }
        }

        // 4. Update State (Expensive Clones happen ONLY here)
        let new_saddle = f64::max(intermediate.saddle_energy, en);
        
        let mut next_remaining = intermediate.remaining_moves.clone();
        next_remaining.remove(cand.original_move_index);

        let mut next_path = intermediate.path.clone();
        next_path.push(cand.applied_move);

        results.push(Intermediate {
            pt: cand.pt,
            saddle_energy: new_saddle,
            current_energy: en,
            remaining_moves: next_remaining,
            path: next_path,
        });
    }

    results
}









// #####################################
//                  TESTS
// #####################################





#[cfg(test)]
mod tests {
    // when importing ff_structure::PairTable into the utils.rs (outside of the tests module), 
    // This warning happens because PairTable is only used inside your tests module.
    // When you run cargo build (or when Rust Analyzer checks your code), 
    // the #[cfg(test)] module is ignored/stripped out. This leaves the import use ff_structure::PairTable; sitting at the top level with no one using it, hence the "unused import" warning.
    use ff_structure::DotBracketVec; 
    use ff_structure::PairTable;
    // use structs defined above (parent module), because tests is a submodule of the module where the structs are defined
    use super::*;

    #[test]
    fn test_valid_pair_table() {
        let pt = PairTable::try_from("((..))").unwrap();
        assert_eq!(pt.len(), 6);
        assert_eq!(pt[0 as usize], Some(5));
        assert_eq!(pt[1 as usize], Some(4));
        assert_eq!(pt[2 as usize], None);
        assert_eq!(pt[3 as usize], None);
        assert_eq!(pt[4 as usize], Some(1));
        assert_eq!(pt[5 as usize], Some(0));
    }

    #[test]
    fn test_compare_structures_and_return_structure_difference() {
        // expected value
        let expected = StructureDifference {
            bp_distance: 7,
            move_list: vec![
                Move { i: 1, j: 20, is_insertion: false }, // del(1, 20)
                Move { i: 2, j: 19, is_insertion: false }, // del(2, 19)
                Move { i: 3, j: 18, is_insertion: false }, // del(3, 18)
                Move { i: 2, j: 18, is_insertion: true },  // ins(2, 18)
                Move { i: 3, j: 17, is_insertion: true },  // ins(3, 17)
                Move { i: 4, j: 16, is_insertion: true },  // ins(4, 16)
                Move { i: 5, j: 15, is_insertion: true },  // ins(5, 15)
            ]};
            

        // call the function under test (replace this with your real call)
        let pt1 = PairTable::try_from(".(((..............)))....").unwrap();
        let pt2 = PairTable::try_from("..((((.........))))......").unwrap();
        println!("pt1: {}", pt1);
        println!("pt2: {}", pt2);
        
        let got = compare_structures(&pt1, &pt2);
        println!("Prepared moves:");
        for m in &got.move_list  {
        println!("{:?}", m);
    }
        //println!("Got StructureDifference: {:?}", got);

        assert_eq!(got, expected);
    }

    #[test]
    fn test_vienna_energy_model(){
        use ff_energy::ViennaRNA;
        // energy_of_structure is exported via the EnergyModel trait
        use ff_energy::EnergyModel;
        use ff_energy::NucleotideVec;


        let model = ViennaRNA::default();

        let seq1 = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let struct1 = ".(((..............)))....";
        let expected_energy = -2.2;
        assert_eq!(model.energy_of_structure(&NucleotideVec::try_from_rna(seq1).expect("failed to read RNA sequence"), &PairTable::try_from(struct1).expect("valid")).expect("failed to calculate energy") as f64/100.0, expected_energy);
        
        let seq = "UCUACUAUUCCGGCUUGACAUAAAUAUCGAGUGCUCGACC";
        let dbr = "...........(.(((((........)))))..)......";
        let exp_energy = -210;
        assert_eq!(model.energy_of_structure(&NucleotideVec::try_from_rna(seq).expect("failed to read RNA sequence"), &PairTable::try_from(dbr).expect("valid")).expect("failed to calculate energy"), exp_energy);
        
    }


    #[test]
    fn test_generate_valid_neighbors_01() {
        let db1 = "..((..))";
        //let db2 = "((......))";

        let pt1 = PairTable::try_from(db1).expect("Invalid structure");
        //let pt2 = PairTable::try_from(db2).expect("Invalid structure");
        let moves = vec![
            Move { i: 0, j: 7, is_insertion: true }, // InValid
            Move { i: 2, j: 7, is_insertion: false }, // Valid
            Move { i: 0, j: 5, is_insertion: true }, // Invalid (crossing)
        ];
        let neighbors = generate_valid_neighbors(&pt1, &moves);
        assert_eq!(neighbors.len(), 1, "Only one valid move should be found");
        assert_eq!(neighbors[0].applied_move.i, 2);
        assert_eq!(neighbors[0].applied_move.j, 7);
        assert_eq!(neighbors[0].applied_move.is_insertion, false);
    }

    #[test]
fn test_generate_valid_neighbors_02() {
    // --- SCENARIO 1: Simple Valid Moves ---
    // Start: "...." (Length 4, all unpaired)
    let pt_empty = PairTable::try_from("....").unwrap();
    
    // Potential moves:
    // 1. Pair (0, 3) -> "(..)"
    // 2. Pair (1, 2) -> ".(.)"
    // Initialize manually
    let moves = vec![
        Move { i: 0, j: 3, is_insertion: true }, 
        Move { i: 1, j: 2, is_insertion: true }, 
    ];

    let neighbors = generate_valid_neighbors(&pt_empty, &moves);

    println!("--- Scenario 1: Empty Structure ---");
    assert_eq!(neighbors.len(), 2, "Both moves should be valid on empty structure");
    // Verify the structures in the neighbors
    // Neighbor 0 should have 0-3 paired
    assert_eq!(neighbors[0].pt[0 as usize], Some(3));

}

    #[test]
    fn test_generate_valid_neighbors_03_pseudoknot(){
        // --- SCENARIO 2: Crossing / Pseudoknot Rejection ---
    // Start: "(...)." (Length 6)
    // Pair exists at 0-4. 
    // We try to pair index 3 (inside the loop) with index 5 (outside the loop).
    let pt_crossing = PairTable::try_from("(...).").unwrap();
    
    // This move would create a crossing link (pseudoknot)
    let bad_moves = vec![Move { i: 3, j: 5, is_insertion: true }];
    
    let neighbors_cross = generate_valid_neighbors(&pt_crossing, &bad_moves);
    
    println!("\n--- Scenario 2: Crossing Check ---");
    println!("Trying to pair (3,5) on structure '(...).'");
    println!("No valid neighbors, as expected. Would create pseudoknot.");
    assert_eq!(neighbors_cross.len(), 0, "Crossing move should be rejected automatically");

    }

    #[test]
    fn test_generate_valid_neighbors_printing_debugging(){
        // Does not test anything, just for visual inspection and debugging of the generated neighbors.
        let db1 = ".(((..............)))....";
        let db2 = "..((((.........))))......";

        let pt1 = PairTable::try_from(db1).expect("Invalid structure");
        let pt2 = PairTable::try_from(db2).expect("Invalid structure");

        let comparison = compare_structures(&pt1, &pt2);
        let move_list = comparison.move_list;

        println!("bp_distance ({}) to go from Structure 1 to Structure 2", move_list.len());
        println!("-----------------------------------------------------");

        let neighbors = generate_valid_neighbors(&pt1, &move_list);
        println!("Generated only {} valid neighbors from Structure 1 towards Structure 2", neighbors.len());
        for (idx, n) in neighbors.iter().enumerate() {
            println!("\nNeighbor {}: {:?} --> New Structure: {}", idx+1, n.applied_move, DotBracketVec::try_from(&n.pt).unwrap());
        }       



    }

    #[test]
    fn test_apply_move() {
        let seq1  = "AGCCAUGAGUGUAUAGUGGGCCUAU";
        let str1 = ".(((..............)))....";
        let str2 = "..((((.........))))......";
        let model = ViennaRNA::default();


        // 1. Setup & Validation
        let seq_vec = NucleotideVec::try_from_rna(seq1).expect("failed to read RNA sequence");
        let pt_start = PairTable::try_from(str1).unwrap();
        let pt_target = PairTable::try_from(str2).unwrap();

        // 2. Prepare Moves
        let comparison = compare_structures(&pt_start, &pt_target);
        let move_list = comparison.move_list;
        

        // 3. Initialize Energies
        let start_energy = model.energy_of_structure(&seq_vec, &pt_start).expect("failed to calculate energy ") as f64 / 100.0;



        // Initialize Intermediate
        let current = Intermediate {
            pt: pt_start,
            saddle_energy: start_energy,
            current_energy: start_energy,
            remaining_moves: move_list,
            path: Vec::new(),
        };
        let expanded = apply_move(&model, &seq_vec, &current, None);
        println!("Expanded into {} new structures.", expanded.len());
        println!("-------------------------------------");
        //println!("Expanded States:");
        for (idx, state) in expanded.iter().enumerate() {
            let db = DotBracketVec::try_from(&state.pt).unwrap();
            println!("State {}: Structure: {}\t Current Energy: {:.2} kcal/mol", idx+1, db, state.current_energy);
        }

        assert!(expanded.len() > 0, "At least one expansion should be possible");


    }



}
