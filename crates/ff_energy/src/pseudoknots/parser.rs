//! Parsing of secondary structure as Extended Dot-Bracket notation to build a PairTable
//! for susequently building the Tree of closed Regions. (closed_region_tree.rs)

use ff_structure::NAIDX;
use ff_structure::{BracketKind, ExtendedDotBracket};
use ff_structure::{PairTable, StructureError};
use std::collections::HashMap;

/// Takes a slice of `ExtendedDotBracket` enum variants and returns the corresponding `PairTable`,
/// or a `StructureError` if brackets are unmatched.
pub fn extended_dot_bracket_to_pair_table(
    edb: &[ExtendedDotBracket],
) -> Result<PairTable, StructureError> {
    // create a new empty PairTable (all None values) of the same length as the input edb
    let mut pt = PairTable::new(edb.len());

    // create a HashMap to hold stacks for each BracketKind to track open brackets
    //NOTE: no need for Option<NAIDX>, since we only push indices of open brackets, and pop them when we encounter a close bracket
    let mut stacks: HashMap<BracketKind, Vec<NAIDX>> = HashMap::new();

    for (i, ch) in edb.iter().enumerate() {
        match ch {
            // Unpaired and break BracketKinds do not affect the stack => continue
            ExtendedDotBracket::Unpaired => {}
            ExtendedDotBracket::Break => {
                return Err(StructureError::InvalidToken(
                    format!("character '{}'", char::from(*ch)),
                    "extended dot-bracket".to_string(),
                    i,
                ));
            }

            // For an open bracket, push the index onto the corresponding stack with the BracketKind as the key
            ExtendedDotBracket::Open(kind) => {
                stacks.entry(*kind).or_default().push(i as NAIDX);
            }
            // For a close bracket, pop from the corresponding stack and record the pairing in the PairTable
            ExtendedDotBracket::Close(kind) => {
                let j = stacks
                    .get_mut(kind) // returns Option<&mut Vec<NAIDX>>
                    .and_then(|s| s.pop()) // if Option type above returns Some(stack), we pass the stack s as a mutable reference and pop from the stack, which returns Option<NAIDX>
                    .ok_or(StructureError::UnmatchedClose(i))?; // if the Option type above returns None, we return an error for unmatched close bracket at index i

                // build the PairTable by setting the indices i and j to point to each other
                pt[i] = Some(j);
                pt[j as usize] = Some(i as NAIDX);
            }
        }
    }

    // Check for any unclosed open brackets across all kinds after the loop.
    // If any stack is not empty, it means there are unmatched open brackets, and we return an error for the first unmatched open bracket found.
    for kind in BracketKind::all() {
        if let Some(stack) = stacks.get(kind)
            && let Some(&i) = stack.first()
        {
            return Err(StructureError::UnmatchedOpen(i as usize));
        }
    }

    Ok(pt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_edb() {
        let edb = vec![];
        let pt = extended_dot_bracket_to_pair_table(&edb).unwrap();
        assert_eq!(pt.len(), 0);
    }

    #[test]
    fn test_unpaired_edb() {
        let edb = vec![
            ExtendedDotBracket::Unpaired,
            ExtendedDotBracket::Unpaired,
            ExtendedDotBracket::Unpaired,
        ];
        let pt = extended_dot_bracket_to_pair_table(&edb).unwrap();
        assert_eq!(pt.len(), 3);
        assert_eq!(pt[0 as NAIDX], None);
        assert_eq!(pt[1 as NAIDX], None);
        assert_eq!(pt[2 as NAIDX], None);
    }

    #[test]
    fn test_invalid_edb_mismatched_brackets() {
        // Example: "(.]" => mismatched brackets
        let edb = vec![
            ExtendedDotBracket::Open(BracketKind::Round),
            ExtendedDotBracket::Unpaired,
            ExtendedDotBracket::Close(BracketKind::Square), // mismatched close bracket
        ];
        let err = extended_dot_bracket_to_pair_table(&edb).unwrap_err();
        assert!(matches!(err, StructureError::UnmatchedClose(2)));
    }

    #[test]
    fn test_invalid_edb_unmatched_open() {
        // Example: "([.)" => unmatched open bracket
        let edb = vec![
            ExtendedDotBracket::Open(BracketKind::Round),
            ExtendedDotBracket::Open(BracketKind::Square), // unmatched open bracket
            ExtendedDotBracket::Unpaired,
            ExtendedDotBracket::Close(BracketKind::Round),
        ];
        let err = extended_dot_bracket_to_pair_table(&edb).unwrap_err();
        // matches! macro checks if the error is of type UnmatchedOpen and that the index of the unmatched open bracket is 1 (the index of the Square open bracket)
        assert!(matches!(err, StructureError::UnmatchedOpen(1)));
    }
    #[test]
    fn test_invalid_edb_character() {
        // Example: "([x])" => invalid character 'x'
        let edb = vec![
            ExtendedDotBracket::Open(BracketKind::Round),
            ExtendedDotBracket::Open(BracketKind::Square),
            ExtendedDotBracket::Break, // invalid character represented Break variant
            ExtendedDotBracket::Close(BracketKind::Square),
            ExtendedDotBracket::Close(BracketKind::Round),
        ];
        let err = extended_dot_bracket_to_pair_table(&edb).unwrap_err();
        // matches! macro checks if the error is of type InvalidToken and wildcard "_" is used to ignore the specific error "token", "src" and "position", since we are only interested in the type of error
        assert!(matches!(err, StructureError::InvalidToken(_, _, _)));
    }

    #[test]
    fn test_valid_pair_table() {
        // Example: "([.])" => 0 pairs with 4, 1 pairs with 3, 2 is unpaired
        let edb = vec![
            ExtendedDotBracket::Open(BracketKind::Round),
            ExtendedDotBracket::Open(BracketKind::Square),
            ExtendedDotBracket::Unpaired,
            ExtendedDotBracket::Close(BracketKind::Square),
            ExtendedDotBracket::Close(BracketKind::Round),
        ];
        let pt = extended_dot_bracket_to_pair_table(&edb).unwrap();
        assert_eq!(pt[0 as NAIDX], Some(4));
        assert_eq!(pt[1 as NAIDX], Some(3));
        assert_eq!(pt[2 as NAIDX], None);
        assert_eq!(pt[3 as NAIDX], Some(1));
        assert_eq!(pt[4 as NAIDX], Some(0));
    }
}
