//! `PairTable` → extended dot-bracket.

use ff_structure::{BracketKind, PairTable, StructureError};

/// Convert a [`PairTable`] to an extended dot-bracket string.
///
/// Bracket families (`()`, `[]`, `{}`, …) are assigned by greedily colouring
/// the pair-crossing graph: crossing pairs receive different families; the
/// lowest available family index is always chosen.
///
/// Returns `Err` if the structure needs more families than `BracketKind` has.
pub fn pair_table_to_dot_bracket(pt: &PairTable) -> Result<String, StructureError> {
    let n = pt.len();

    // Collect canonical (i < j) pairs sorted by left endpoint.
    let mut pairs: Vec<(usize, usize)> = (0..n)
        .filter_map(|i| pt[i].map(|j| (i, j as usize)))
        .filter(|&(i, j)| i < j)
        .collect();
    pairs.sort_unstable_by_key(|&(i, _)| i);

    let max_levels = BracketKind::all().len();
    let mut family = vec![0u8; pairs.len()];
    let mut blocked = vec![false; max_levels]; // one slot per BracketKind

    for idx in 0..pairs.len() {
        let (pi, pj) = pairs[idx];
        blocked.fill(false);

        // All previously assigned pairs have left endpoint qi <= pi (sort order).
        // (pi, pj) and (qi, qj) cross iff qi < pi < qj < pj, i.e. pi < qj < pj.
        for jdx in 0..idx {
            let (_, qj) = pairs[jdx];
            if pi < qj && qj < pj {
                blocked[family[jdx] as usize] = true;
            }
        }

        let f = (0..max_levels).find(|&f| !blocked[f]).ok_or_else(|| {
            StructureError::InvalidToken(
                format!("pair ({pi},{pj})"),
                "extended dot-bracket".to_string(),
                pi,
            )
        })?;
        family[idx] = f as u8;
    }

    let kinds = BracketKind::all();
    let mut chars: Vec<char> = vec!['.'; n];
    for (idx, &(i, j)) in pairs.iter().enumerate() {
        let kind = kinds[family[idx] as usize];
        chars[i] = kind.open_char();
        chars[j] = kind.close_char();
    }
    Ok(chars.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pseudoknots::extended_dot_bracket_to_pair_table;
    use ff_structure::ExtendedDotBracketVec;

    #[test]
    fn roundtrip_simple_hairpin() {
        let pt = PairTable::try_from("(((...)))").unwrap();
        assert_eq!(pair_table_to_dot_bracket(&pt).unwrap(), "(((...)))");
    }

    #[test]
    fn roundtrip_h_type_pk() {
        let edb = ExtendedDotBracketVec::try_from("(((...[[[...)))...]]]").unwrap();
        let pt = extended_dot_bracket_to_pair_table(&edb).unwrap();
        assert_eq!(
            pair_table_to_dot_bracket(&pt).unwrap(),
            "(((...[[[...)))...]]]"
        );
    }

    #[test]
    fn roundtrip_fully_unpaired() {
        let pt = PairTable::try_from(".....").unwrap();
        assert_eq!(pair_table_to_dot_bracket(&pt).unwrap(), ".....");
    }

    /// `n` single-pair bands, all mutually crossing: open in family order, then
    /// close in the same order.
    fn crossing(n: usize) -> String {
        let kinds = BracketKind::all();
        (0..n)
            .map(|k| kinds[k % kinds.len()].open_char())
            .chain((0..n).map(|k| kinds[k % kinds.len()].close_char()))
            .collect()
    }

    /// Every family, which also pins the greedy colouring to `all()`'s order.
    #[test]
    fn roundtrip_every_family() {
        let s = crossing(BracketKind::all().len());
        let edb = ExtendedDotBracketVec::try_from(s.as_str()).unwrap();
        let pt = extended_dot_bracket_to_pair_table(&edb).unwrap();
        assert_eq!(pair_table_to_dot_bracket(&pt).unwrap(), s);
    }

    /// One band past the alphabet. Built by hand since it has no notation.
    #[test]
    fn one_family_past_the_alphabet() {
        let n = BracketKind::all().len() + 1;
        let mut pt = PairTable::new(2 * n);
        for i in 0..n {
            pt[i] = Some((i + n) as u16);
            pt[i + n] = Some(i as u16);
        }
        assert!(pair_table_to_dot_bracket(&pt).is_err());
    }
}
