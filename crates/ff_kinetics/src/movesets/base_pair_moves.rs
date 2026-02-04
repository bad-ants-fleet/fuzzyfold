use ff_structure::NAIDX;

#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
pub enum Move {
    /// A pair formation.
    Add {
        i: NAIDX,
        j: NAIDX,
    },
    /// A pair opening.
    Del {
        i: NAIDX,
        j: NAIDX,
    },
    /// A three-way shift move (i, j) -> (i, k)
    ShiftIK {
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
    },
    /// A three-way shift move (i, j) -> (j, k)
    ShiftJK {
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
    },
    /// A four-way shift move with adjacent pairs: 
    /// (i, j), (k, l) -> (i, l), (j, k)
    ShiftILJK {
        // Outer base-pair
        i: NAIDX,
        j: NAIDX,
        // Inner base-pair
        k: NAIDX,
        l: NAIDX,
    },
    /// A four-way shift move with enclosing pair: 
    /// (i, j), (k, l) -> (i, k), (l, j)
    ShiftIKLJ {
        // Outer base-pair
        i: NAIDX,
        j: NAIDX,
        // Inner base-pair
        k: NAIDX,
        l: NAIDX,
    },
}

impl Move {
    /// Returns the inverse move.
    pub fn inverse(self) -> Self {
        match self {
            Move::Add { i, j } => Move::Del { i, j },
            Move::Del { i, j } => Move::Add { i, j },
            Move::ShiftIK { i, j, k } => {
                if i < k {
                    Move::ShiftIK { i, j: k, k: j }
                } else {
                    Move::ShiftJK { i: k, j: i, k: j }
                }
            }
            Move::ShiftJK { i, j, k } => {
                if k < j {
                    Move::ShiftJK { i: k, j, k: i }
                } else {
                    Move::ShiftIK { i: j, j: k, k: i }
                }
            }
            Move::ShiftILJK { i, j, k, l } => {
                Move::ShiftIKLJ { i, j: l, k: j, l: k }
            }
            Move::ShiftIKLJ { i, j, k, l } => {
                Move::ShiftILJK { i, j: k, k: l, l:j }
            }
        }
    }

    /// Returns the pair added by the move.
    pub fn added_pair(&self) -> (NAIDX, NAIDX) {
        match &self {
            Move::Add { i, j } =>  (*i, *j),
            Move::ShiftIK { i, k, .. } => if i < k { (*i, *k) } else { (*k, *i) },
            Move::ShiftJK { j, k, .. } => if j < k { (*j, *k) } else { (*k, *j) },
            _ => unreachable!("Move should add a single pair!"),
        }
    }

    /// Returns the two pairs added by the move.
    pub fn added_pairs(&self) -> ((NAIDX, NAIDX), (NAIDX, NAIDX)) {
        match &self {
            Move::ShiftIKLJ { i, j, k, l } => ((*i, *k), (*l, *j)),
            Move::ShiftILJK { i, j, k, l } => ((*i, *l), (*j, *k)),
            _ => unreachable!("Move should add exactly two pairs!"),
        }
    }

    /// Returns the pair deleted by the move.
    pub fn deleted_pair(&self) -> (NAIDX, NAIDX) {
        match &self {
            Move::Del { i, j } =>  (*i, *j),
            Move::ShiftIK { i, j, .. } => (*i, *j),
            Move::ShiftJK { i, j, .. } => (*i, *j),
            _ => unreachable!("Move should delete a single pair!"),
        }
    }

    /// Returns the two pairs deleted by the move.
    pub fn deleted_pairs(&self) -> ((NAIDX, NAIDX), (NAIDX, NAIDX)) {
        match &self {
            Move::ShiftIKLJ { i, j, k, l } => ((*i, *j), (*k, *l)),
            Move::ShiftILJK { i, j, k, l } => ((*i, *j), (*k, *l)),
            _ => unreachable!("Move should delete exactly two pairs!"),
        }
    }

    /// Reports if there is an overlap between the added pairs in a move
    /// and the given pair.
    pub fn conflicts(&self, pair: (NAIDX, NAIDX)) -> bool {
        #[inline]
        fn overlaps(a: (NAIDX, NAIDX), b: (NAIDX, NAIDX)) -> bool {
            let (p, q) = a;
            let (m, n) = b;
            !(q < m || n < p || (p < m && n < q) || (m < p && q < n))
        }
        match &self {
            Move::Add { .. } => overlaps(self.added_pair(), pair),
            Move::ShiftIK { .. } => overlaps(self.added_pair(), pair),
            Move::ShiftJK { .. } => overlaps(self.added_pair(), pair),
            Move::ShiftILJK { .. } => {
                // TODO: maybe check one and debug_assert the other?
                let (p1, p2) = self.added_pairs();
                overlaps(p1, pair) || overlaps(p2, pair)
            }
            Move::ShiftIKLJ { .. } => {
                // TODO: maybe check one and debug_assert the other?
                let (p1, p2) = self.added_pairs();
                overlaps(p1, pair) || overlaps(p2, pair)
            }
            _ => unreachable!("conflicts() is only valid for moves that add pairs"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pair(a: NAIDX, b: NAIDX) -> (NAIDX, NAIDX) {
        if a < b { (a, b) } else { (b, a) }
    }

    #[test]
    fn test_inverse() {
        let moves = [
            Move::Add { i: 2, j: 10 },
            Move::Del { i: 3, j: 7 },
            Move::ShiftIK { i: 2, j: 8, k: 5 },
            Move::ShiftJK { i: 2, j: 8, k: 5 },
            Move::ShiftIKLJ { i: 1, j: 10, k: 4, l: 7 },
            Move::ShiftILJK { i: 1, j: 10, k: 4, l: 7 },
        ];

        for mv in moves {
            assert_eq!(mv.inverse().inverse(), mv);
        }
    }

    #[test]
    fn add_and_del_pairs() {
        let add = Move::Add { i: 3, j: 9 };
        assert_eq!(add.added_pair(), pair(3, 9));

        let del = Move::Del { i: 3, j: 9 };
        assert_eq!(del.deleted_pair(), (3, 9));
    }

    #[test]
    fn three_way_shift_pairs() {
        let mv = Move::ShiftIK { i: 2, j: 8, k: 5 };
        assert_eq!(mv.deleted_pair(), (2, 8));
        assert_eq!(mv.added_pair(), pair(2, 5));

        let mv = Move::ShiftJK { i: 2, j: 8, k: 5 };
        assert_eq!(mv.deleted_pair(), (2, 8));
        assert_eq!(mv.added_pair(), pair(8, 5));
    }

    #[test]
    fn four_way_shift_iklj_pairs() {
        let mv = Move::ShiftIKLJ { i: 1, j: 10, k: 4, l: 7 };

        let (d1, d2) = mv.deleted_pairs();
        assert_eq!(pair(d1.0, d1.1), pair(1, 10));
        assert_eq!(pair(d2.0, d2.1), pair(4, 7));

        let (a1, a2) = mv.added_pairs();
        assert_eq!(pair(a1.0, a1.1), pair(1, 4));
        assert_eq!(pair(a2.0, a2.1), pair(7, 10));
    }

    #[test]
    fn four_way_shift_iljk_pairs() {
        let mv = Move::ShiftILJK { i: 1, j: 10, k: 4, l: 7 };

        let (d1, d2) = mv.deleted_pairs();
        assert_eq!(pair(d1.0, d1.1), pair(1, 10));
        assert_eq!(pair(d2.0, d2.1), pair(4, 7));

        let (a1, a2) = mv.added_pairs();
        assert_eq!(pair(a1.0, a1.1), pair(1, 7));
        assert_eq!(pair(a2.0, a2.1), pair(4, 10));
    }

}
