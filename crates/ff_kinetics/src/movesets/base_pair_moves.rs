use ff_structure::NAIDX;

#[derive(Debug, Clone, Copy, Hash, Eq, PartialEq)]
pub enum Move {
    Add {
        i: NAIDX,
        j: NAIDX,
    },
    Del {
        i: NAIDX,
        j: NAIDX,
    },
    ShiftIK {
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
    },
    ShiftJK {
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
    },
    //FourwayShift {
    //    i: NAIDX,
    //    j: NAIDX,
    //    k: NAIDX,
    //    l: NAIDX,
    //},
}

impl Move {
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
        }
    }

    pub fn added_pair(&self) -> (NAIDX, NAIDX) {
        match &self {
            Move::Add { i, j } =>  (*i, *j),
            Move::ShiftIK { i, k, .. } => if i < k { (*i, *k) } else { (*k, *i) },
            Move::ShiftJK { j, k, .. } => if j < k { (*j, *k) } else { (*k, *j) },
            _ => unreachable!(""),
        }
    }

    pub fn deleted_pair(&self) -> (NAIDX, NAIDX) {
        match &self {
            Move::Del { i, j } =>  (*i, *j),
            Move::ShiftIK { i, j, .. } => (*i, *j),
            Move::ShiftJK { i, j, .. } => (*i, *j),
            _ => unreachable!(""),
        }
    }

    pub fn conflicts(&self, pair: (NAIDX, NAIDX)) -> bool {

        #[inline]
        fn overlaps(a: (NAIDX, NAIDX), b: (NAIDX, NAIDX)) -> bool {
            let (p, q) = a;
            let (m, n) = b;
            !(q < m || n < p || (p < m && n < q) || (m < p && q < n))
        }

        overlaps(self.added_pair(), pair)
    }
}


