/// The shift moves a move set proposes, besides adding and deleting pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ShiftMoves {
    /// A pair moves one of its bases: (i, j) -> (i, k) or (j, k).
    pub three_way: bool,
    /// Two pairs exchange their bases: (i, j), (k, l) -> (i, l), (j, k) or
    /// (i, k), (l, j).
    pub four_way: bool,
}

impl ShiftMoves {
    /// Only adding and deleting pairs.
    pub const NONE: ShiftMoves = ShiftMoves { three_way: false, four_way: false };
    /// Three-way and four-way shifts.
    pub const ALL: ShiftMoves = ShiftMoves { three_way: true, four_way: true };
}
