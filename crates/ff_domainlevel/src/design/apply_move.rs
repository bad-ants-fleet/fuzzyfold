
use ff_structure::Pair;
use ff_structure::PairTable;
use ff_structure::LoopTable;
use ff_structure::LoopInfo;

pub trait ApplyMove {
    fn try_move(&self, pair: Pair) -> Result<Option<Pair>, String>;
    fn apply_move(&mut self, old: Option<Pair>, new: Pair);
}

impl ApplyMove for PairTable {

    fn try_move(&self, pair: Pair) -> Result<Option<Pair>, String> {
        use LoopInfo::*;
        let (i, j) = (pair.i() as usize, pair.j() as usize);
        if Some(pair.j()) == self[i] && Some(pair.i()) == self[j] {
            return Ok(Some(pair));
        }
        let lt = LoopTable::from(self);
        match (lt[i], lt[j]) {
            (Unpaired { l: iloop }, Unpaired { l: jloop }) => {
                if iloop == jloop {
                    Ok(None)
                } else {
                    Err("Unpaired bases are in different loops.".to_string())
                }
            }
            (Unpaired { l: iloop }, Paired { i: inner_loop, o: outer_loop }) => {
                if iloop == inner_loop || iloop == outer_loop {
                    let pi = self[j].unwrap();
                    if pi < pair.j() { 
                        Ok(Some(Pair::new(pi, pair.j())))
                    } else { 
                        Ok(Some(Pair::new(pair.j(), pi)))
                    }
                } else {
                    Err(format!("Loop mismatch ({i} unpaired, {j} paired)."))
                }
            }
            (Paired { i: inner_loop, o: outer_loop }, Unpaired { l: jloop }) => {
                if jloop == inner_loop || jloop == outer_loop {
                    let pj = self[i].unwrap();
                    if pj < pair.i() { 
                        Ok(Some(Pair::new(pj, pair.i())))
                    } else { 
                        Ok(Some(Pair::new(pair.i(), pj)))
                    }
                } else {
                    Err(format!("Loop mismatch ({i} paired, {j} unpaired)."))
                }
            }
            (Paired { i: i_inner, o: i_outer }, Paired { i: j_inner, o: j_outer }) => {
                debug_assert!(i_inner != j_inner);
                if i_outer == j_outer || i_outer == j_inner || j_outer == i_inner {
                    Err(format!("Both bases paired, but could work. {i} {j}"))
                } else {
                    Err(format!("Both bases paired and loop mismatch! {}", pair))
                }
            }
        }
    }

    fn apply_move(&mut self, old: Option<Pair>, new: Pair) {
        if let Some(old) = old {
            self[old.i() as usize] = None;
            self[old.j() as usize] = None;
        }
        self[new.i() as usize] = Some(new.j());
        self[new.j() as usize] = Some(new.i());
    }

}

