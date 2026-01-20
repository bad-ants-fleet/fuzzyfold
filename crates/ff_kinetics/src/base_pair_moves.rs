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
    ShiftJ {
        i: NAIDX,
        j: NAIDX,
        k: NAIDX,
    },
    ShiftI {
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
            Move::ShiftJ { i, j, k } => {
                if i < k {
                    Move::ShiftJ { i, k, j }
                } else {
                    Move::ShiftI { k, i, j }
                }
            }
            Move::ShiftI { i, j, k } => {
                if k < j {
                    Move::ShiftI { k, j, i }
                } else {
                    Move::ShiftJ { j, k, i }
                }
           }
        }
    }
}


