//! MultiPairTable: the pairs of a structure of one or more strands.

use std::fmt;
use std::convert::TryFrom;

use crate::error::StructureError;
use crate::dotbracket::DotBracket;
use crate::dotbracket::DotBracketVec;

/// How a position is stored in the table. `u16` keeps the table small and
/// limits a structure to `NAIDX::MAX - 2` positions; change it here for
/// longer molecules.
#[allow(clippy::upper_case_acronyms)]
type NAIDX = u16;

/// The state of one position of a [`MultiPairTable`]: paired with the
/// given position, unpaired, or the end of a strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MultiStruct {
    Paired(usize),
    Unpaired,
    /// `+`: the strand before it ends.
    StrandBreak,
    /// `^`: the strand before it is circular, it continues at its own first
    /// base.
    StrandWrap,
}

/// Stored in a position that is unpaired.
const UNPAIRED: NAIDX = NAIDX::MAX;
/// Stored in a position that is a strand break.
const STRAND_BREAK: NAIDX = NAIDX::MAX - 1;
/// Stored in a position that closes a circular strand. The smallest reserved
/// value: positions are below it.
const STRAND_WRAP: NAIDX = NAIDX::MAX - 2;

/// The pairs of a structure of one or more strands, with or without
/// pseudoknots. This is the input of the loop decomposition and energy
/// evaluation in `ff_energy` (`RingDecomposition`, `EnergyModel`), of the
/// kinetics in `ff_kinetics`, and of [`ClosedRegions`](crate::ClosedRegions).
///
/// Each position stores its partner, or a value reserved for unpaired 
/// positions and the ends of strands, in two bytes. A structure has at most
/// 65533 positions, strand ends included.
///
/// A strand ends at a strand break (`+`), at a strand wrap (`^`), or at
/// the end of the structure. A strand wrap joins the strand before it to
/// its own first base, so `(((...)))^` is a circular molecule. In the order
/// of positions, a strand wrap separates strands like a strand break.
///
/// # Parse and write
///
/// ```
/// use ff_structure::{DotBracketVec, MultiPairTable};
///
/// let mpt = MultiPairTable::try_from("((.[[+))..]]").unwrap();
/// let db = DotBracketVec::try_from(&mpt).unwrap();
/// assert_eq!(db.to_string(), "((.[[+))..]]");
/// ```
///
/// # Read positions
///
/// [`get`](Self::get) returns the [`MultiStruct`] of a position;
/// [`iter`](Self::iter) visits all positions 5' to 3'.
///
/// ```
/// use ff_structure::{MultiPairTable, MultiStruct};
///
/// let mpt = MultiPairTable::try_from("((.[[+))..]]").unwrap();
/// assert_eq!(mpt.get(0), MultiStruct::Paired(7));
/// assert_eq!(mpt.get(2), MultiStruct::Unpaired);
///
/// // The pairs (i, j) with i < j.
/// let pairs: Vec<(usize, usize)> = mpt.iter().enumerate()
///     .filter_map(|(i, s)| match s {
///         MultiStruct::Paired(j) if i < j => Some((i, j)),
///         _ => None,
///     })
///     .collect();
/// assert_eq!(pairs, vec![(0, 7), (1, 6), (3, 11), (4, 10)]);
///
/// // The strand breaks.
/// let breaks: Vec<usize> = (0..mpt.len())
///     .filter(|&i| mpt.get(i) == MultiStruct::StrandBreak)
///     .collect();
/// assert_eq!(breaks, vec![5]);
/// ```
///
/// # Build from pairs
///
/// [`set`](Self::set) changes one position and does not check the partner,
/// so set both ends of a pair.
///
/// ```
/// use ff_structure::{DotBracketVec, MultiPairTable, MultiStruct};
///
/// let mut mpt = MultiPairTable::try_from("....+....").unwrap();
/// for (i, j) in [(0, 8), (1, 7)] {
///     mpt.set(i, MultiStruct::Paired(j));
///     mpt.set(j, MultiStruct::Paired(i));
/// }
/// assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), "((..+..))");
/// ```
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct MultiPairTable(Vec<NAIDX>);

impl MultiPairTable {
    /// The state of position `i`. Panics if `i` is out of bounds.
    pub fn get(&self, i: usize) -> MultiStruct {
        match self.0[i] {
            UNPAIRED => MultiStruct::Unpaired,
            STRAND_BREAK => MultiStruct::StrandBreak,
            STRAND_WRAP => MultiStruct::StrandWrap,
            p => MultiStruct::Paired(p as usize),
        }
    }

    /// Set the state of position `i`. The partner is not updated. Panics if
    /// `i` is out of bounds, or the partner is too large for the table.
    pub fn set(&mut self, i: usize, state: MultiStruct) {
        self.0[i] = match state {
            MultiStruct::Unpaired => UNPAIRED,
            MultiStruct::StrandBreak => STRAND_BREAK,
            MultiStruct::StrandWrap => STRAND_WRAP,
            MultiStruct::Paired(p) => {
                assert!(p < STRAND_WRAP as usize, "Position {} is out of range.", p);
                p as NAIDX
            }
        };
    }

    /// True if no two pairs cross: a structure without pseudoknots. A strand
    /// wrap counts as a pair of the first and the last base of its circle,
    /// so a pair from a circle to another strand crosses.
    pub fn is_nested(&self) -> bool {
        let mut open = Vec::new();
        // The depth at the start of the strand, and the lowest depth on it.
        let (mut start, mut low) = (0, 0);
        for k in 0..self.len() {
            match self.get(k) {
                MultiStruct::Paired(p) if k < p => open.push(k),
                MultiStruct::Paired(p) => {
                    if open.pop() != Some(p) {
                        return false;
                    }
                    low = low.min(open.len());
                }
                MultiStruct::StrandWrap if low < start || open.len() != start => return false,
                MultiStruct::StrandBreak | MultiStruct::StrandWrap => {
                    start = open.len();
                    low = start;
                }
                MultiStruct::Unpaired => (),
            }
        }
        true
    }

    /// The number of positions, strand ends included.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True for a table without positions.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// The state of each position, 5' to 3'.
    pub fn iter(&self) -> impl Iterator<Item = MultiStruct> + '_ {
        (0..self.len()).map(|i| self.get(i))
    }
}

impl fmt::Debug for MultiPairTable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

/// Brackets of different levels match independently. An unmatched opening
/// bracket is reported at the first such position. Every strand has at
/// least one position: a strand break or strand wrap must follow one,
/// and a structure must not end with a strand break.
impl TryFrom<&DotBracketVec> for MultiPairTable {
    type Error = StructureError;

    fn try_from(db: &DotBracketVec) -> Result<Self, Self::Error> {
        if db.len() > STRAND_WRAP as usize {
            return Err(StructureError::TooLong(db.len(), STRAND_WRAP as usize));
        }
        // Per level: the opening brackets that are not matched yet.
        let mut open: Vec<Vec<usize>> = vec![Vec::new(); DotBracket::LEVELS];
        let mut table = vec![UNPAIRED; db.len()];
        for (j, &d) in db.iter().enumerate() {
            if d == DotBracket::STRAND_BREAK || d == DotBracket::STRAND_WRAP {
                let empty = j == 0 || matches!(table[j - 1], STRAND_BREAK | STRAND_WRAP);
                if empty || (d == DotBracket::STRAND_BREAK && j + 1 == db.len()) {
                    return Err(StructureError::InvalidToken(
                        "empty strand".into(), "a structure".into(), j));
                }
                table[j] = if d == DotBracket::STRAND_BREAK { STRAND_BREAK } else { STRAND_WRAP };
            } else if let Some(level) = d.open_level() {
                open[level].push(j);
            } else if let Some(level) = d.close_level() {
                let i = open[level].pop()
                    .ok_or(StructureError::UnmatchedClose(d.into(), j))?;
                table[i] = j as NAIDX;
                table[j] = i as NAIDX;
            }
        }
        if let Some(&i) = open.iter().filter_map(|stack| stack.first()).min() {
            return Err(StructureError::UnmatchedOpen(db[i].into(), i));
        }
        Ok(MultiPairTable(table))
    }
}

impl TryFrom<&str> for MultiPairTable {
    type Error = StructureError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        MultiPairTable::try_from(&DotBracketVec::try_from(s)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use MultiStruct::{StrandWrap, Paired, StrandBreak, Unpaired};

    #[test]
    fn nested_with_circles() {
        let nested = |db: &str| MultiPairTable::try_from(db).unwrap().is_nested();
        for db in ["(((...)))", "((.+))", "((...))^", "((...))^((...))", "((..+....^..))", "(..+.(.).^..)", "..^..^"] {
            assert!(nested(db), "{db}");
        }
        for db in ["((.[[.))..]]", "((..^..))", "(..^..)", "(..+..(^..))", "..(..^..)", "(..^)..+.."] {
            assert!(!nested(db), "{db}");
        }
    }

    #[test]
    fn test_valid_pair_table() {
        let mpt = MultiPairTable::try_from("((.+))").unwrap();
        assert_eq!(mpt.len(), 6);
        assert_eq!(mpt.iter().collect::<Vec<_>>(),
            vec![Paired(5), Paired(4), Unpaired, StrandBreak, Paired(1), Paired(0)]);
    }

    #[test]
    fn test_strand_wrap() {
        // A circle bound by a linear strand.
        let mpt = MultiPairTable::try_from("((((....^))))").unwrap();
        assert_eq!(mpt.get(8), StrandWrap);
        assert_eq!(mpt.get(0), Paired(12));
        assert_eq!(DotBracketVec::try_from(&mpt).unwrap().to_string(), "((((....^))))");
        for s in ["^(..)", "(..)+^", "(..)^^", "+(..)", "(.++.)", "(..)+"] {
            let err = MultiPairTable::try_from(s).unwrap_err();
            assert!(err.to_string().starts_with("Invalid empty strand"), "{s}");
        }
    }

    #[test]
    fn test_pseudoknot() {
        let mpt = MultiPairTable::try_from("([)]").unwrap();
        assert_eq!(mpt.iter().collect::<Vec<_>>(),
            vec![Paired(2), Paired(3), Paired(0), Paired(1)]);
    }

    #[test]
    fn test_parse_errors() {
        for (s, msg) in [
            ("(#)", "Invalid character '#' in dot-bracket at position 1"),
            ("(()", "Unmatched '(' at position 0"),
            ("())", "Unmatched ')' at position 2"),
            ("[[]", "Unmatched '[' at position 0"),
            ("[]]", "Unmatched ']' at position 2"),
            // Both are unmatched: the first one is reported, whatever its level.
            ("[(", "Unmatched '[' at position 0"),
        ] {
            assert_eq!(MultiPairTable::try_from(s).unwrap_err().to_string(), msg, "{s}");
        }
    }

    #[test]
    fn test_too_long() {
        let s = ".".repeat(STRAND_WRAP as usize + 1);
        let err = MultiPairTable::try_from(s.as_str()).unwrap_err();
        assert!(matches!(err, StructureError::TooLong(n, max) if n == max + 1 && max == STRAND_WRAP as usize));
        assert!(MultiPairTable::try_from(&s[1..]).is_ok());
    }

}
