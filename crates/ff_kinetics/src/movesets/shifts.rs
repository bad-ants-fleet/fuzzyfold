//! Shift moves: three-way shifts (a pair moves one of its bases) and
//! four-way shifts (two pairs exchange their bases), with their barriers.


use ff_energy::NestedLoop;
use ff_energy::EnergyModel;

use crate::Move;
use crate::Neighbors;
use crate::LoopTable;

type Pair = (usize, usize);


/// The two loops on either side of a pair, (outer, inner), and the loop
/// without the pair: a three-way shift of the pair splits it again.
struct Joined {
    pair: Pair,
    outer: usize,
    inner: usize,
    combo: NestedLoop,
}

impl Joined {
    fn of<E: EnergyModel>(ltab: &LoopTable<E>, (i, j): Pair) -> Self {
        let (outer, inner) = (ltab.loop_lookup(i), ltab.loop_lookup(j));
        let combo = ltab.get(outer).0.join_loop(&ltab.get(inner).0);
        Joined { pair: (i, j), outer, inner, combo }
    }
}

/// The three-way shift moves into the unpaired bases of a loop, with
/// their barriers.
///
/// A three-way shift turns the two loops on either side of a pair (i, j)
/// into two new loops: (outer0, inner0) -> (outer1, inner1). Its barrier is
/// the largest of
///  - (E(outer1) + E(inner1)) - (E(outer0) + E(inner0)),
///  - E(outer1) - E(outer0),
///  - E(inner1) - E(inner0).
pub(crate) fn three_way_shifts<E: EnergyModel>(ltab: &LoopTable<E>, index: usize) -> Neighbors {
    let (init, _) = ltab.get(index);
    // Each pair joined out once, for all of its shifts.
    let joined: Vec<Joined> = init.index_pairs().into_iter()
        .map(|pair| Joined::of(ltab, pair))
        .collect();

    let mut neighbors: Neighbors = Vec::new();
    for s in init.stretches() {
        let (p5, p3) = s.ends;
        for k in s.unpaired {
            for jn in &joined {
                let (i, j) = jn.pair;
                for p in [i, j] {
                    if !ltab.can_pair(k, p) {
                        continue;
                    }
                    if (p == p5 && p + ltab.min_hairpin_size() >= k) 
                        || (p == p3 && k + ltab.min_hairpin_size() >= p) 
                    {
                        continue;
                    } 
                    let mv = if p == i {
                        Move::ShiftIK { i, j, k }
                    } else {
                        Move::ShiftJK { i, j, k }
                    };
                    if let Some(delta) = three_way_barrier(ltab, jn, &mv) {
                        neighbors.push((mv, delta))
                    }
                }
            }
        }
    }
    neighbors
}
 
/// The barrier of a three-way shift: the largest of the energy change
/// and the changes of either loop. None if a new loop does not fit (see
/// `LoopTable::fits`).
fn three_way_barrier<E: EnergyModel>(
    ltab: &LoopTable<E>,
    joined: &Joined,
    mv: &Move,
) -> Option<i32> {
    debug_assert_eq!(joined.pair, mv.deleted_pair());
    let (outer0_idx, inner0_idx) = (joined.outer, joined.inner);
    let (p, q) = mv.added_pair();
    let (outer1, inner1) = joined.combo.split_loop(p, q);
    if !ltab.fits(&outer1) || !ltab.fits(&inner1) {
        return None;
    }

    let (_o, outer0_en) = ltab.get(outer0_idx);
    let (_i, inner0_en) = ltab.get(inner0_idx);
    let outer1_en = ltab.energy_of_loop(&outer1);
    let inner1_en = ltab.energy_of_loop(&inner1);

    #[cfg(feature = "shift_analysis")]
    print_three_way(_o, _i, &outer1, &inner1, (*outer0_en, *inner0_en, outer1_en, inner1_en));

    Some(((outer1_en + inner1_en) - (outer0_en + inner0_en))
        .max(outer1_en - outer0_en)
        .max(inner1_en - inner0_en))
}

/// The four-way shift moves of the pairs of a loop, with their barriers.
///
/// A four-way shift turns three loops into three new loops: (center,
/// merge1, merge2) -> (inside, outer1, outer2), where inside joins the two
/// merged loops and the center splits into the two outer loops. Its
/// barrier is the largest of
///  - (E(inside) + E(outer1) + E(outer2)) - (E(center) + E(merge1) + E(merge2)),
///  - E(outer1) + E(outer2) - E(center),
///  - E(inside) - E(merge1) - E(merge2).
pub(crate) fn four_way_shifts<E: EnergyModel>(ltab: &LoopTable<E>, index: usize) -> Neighbors {
    let (init, _) = ltab.get(index);
    let pairs = init.index_pairs();

    let mut neighbors = Vec::new(); 
    for (p1, &(i, j)) in pairs.iter().enumerate() {
        for &(k, l) in pairs[p1 + 1..].iter() {
            if j + ltab.min_hairpin_size() < k 
                && ltab.can_pair(i, l) 
                && ltab.can_pair(j, k) 
            {
                let mv = Move::ShiftILJK { i, j, k, l };
                if let Some(delta) = four_way_barrier(ltab, &mv) {
                    neighbors.push((mv, delta));
                }
            } else if ltab.can_pair(i, k) && ltab.can_pair(l, j) 
                && i + ltab.min_hairpin_size() < k 
                && l + ltab.min_hairpin_size() < j
            {
                let mv = Move::ShiftIKLJ { i, j, k, l };
                if let Some(delta) = four_way_barrier(ltab, &mv) {
                    neighbors.push((mv, delta));
                }
            }
        }
    }
    neighbors
}

/// The loops of a four-way shift: the indices of the three loops it
/// changes, and the three new loops.
pub(crate) fn four_way_loops<E: EnergyModel>(
    ltab: &LoopTable<E>,
    mv: &Move,
) -> (usize, // center
      usize, // merge1
      usize, // merge2
      NestedLoop, // inside
      NestedLoop, // outer1
      NestedLoop) // outer2
{
    match mv {
        Move::ShiftIKLJ { i, j, k, l } => {
            let merge1_idx = ltab.loop_lookup(*i);
            let center_idx = ltab.loop_lookup(*j);
            let merge2_idx = ltab.loop_lookup(*l);
            let (merge1, _energy) = ltab.get(merge1_idx);
            let (center, center_en) = ltab.get(center_idx);
            debug_assert_eq!(center_en, &ltab.geti(*k).1);
            let (merge2, _energy) = ltab.get(merge2_idx);

            let inside = merge1.join_loop(center).join_loop(merge2);
            let (inside, outer1) = inside.split_loop(*i, *k);
            let (inside, outer2) = inside.split_loop(*l, *j);
            (center_idx, merge1_idx, merge2_idx, inside, outer1, outer2)
        },
        Move::ShiftILJK { i, j, k, l } => {

            let center_idx = ltab.loop_lookup(*i);
            let merge1_idx = ltab.loop_lookup(*j);
            let merge2_idx = ltab.loop_lookup(*l);
            let (merge1, _energy) = ltab.get(merge1_idx);
            let (center, center_en) = ltab.get(center_idx);
            debug_assert_eq!(center_en, &ltab.geti(*k).1);
            let (merge2, _energy) = ltab.get(merge2_idx);

            let inside = center.join_loop(merge1).join_loop(merge2);
            let (outer1, inside) = inside.split_loop(*i, *l);
            let (inside, outer2) = inside.split_loop(*j, *k);
            (center_idx, merge1_idx, merge2_idx, inside, outer1, outer2)
        },
        _ => panic!("FourWayNeighbors called with wrong move {:?}", mv),
    }
}
 
/// The barrier of a four-way shift: the largest of the energy change, the
/// change of the center loop split in two, and of the merged loops. None
/// if a new loop does not fit (see `LoopTable::fits`).
fn four_way_barrier<E: EnergyModel>(
    ltab: &LoopTable<E>,
    mv: &Move,
) -> Option<i32> {
    let (center_idx, merge1_idx, merge2_idx, 
         inside, outer1, outer2) = four_way_loops(ltab, mv);
    if !ltab.fits(&inside) || !ltab.fits(&outer1) || !ltab.fits(&outer2) {
        return None;
    }

    let (_, center_en) = ltab.get(center_idx);
    let (_m1, merge1_en) = ltab.get(merge1_idx);
    let (_m2, merge2_en) = ltab.get(merge2_idx);

    let inside_en = ltab.energy_of_loop(&inside);
    let outer1_en = ltab.energy_of_loop(&outer1);
    let outer2_en = ltab.energy_of_loop(&outer2);

    #[cfg(feature = "shift_analysis")]
    print_four_way(_m1, _m2, &outer1, &outer2,
        (*center_en, *merge1_en, *merge2_en, inside_en, outer1_en, outer2_en));

    Some(((inside_en + outer1_en + outer2_en) - (center_en + merge1_en + merge2_en))
        .max(outer1_en + outer2_en - center_en)
        .max(inside_en - merge1_en - merge2_en))
}

/// For a figure in Badelt et al. 2026: the energies of a three-way shift and
/// the kind of rearrangement, printed.
#[cfg(feature = "shift_analysis")]
fn print_three_way(_o: &NestedLoop, _i: &NestedLoop, outer1: &NestedLoop, inner1: &NestedLoop,
    (outer0_en, inner0_en, outer1_en, inner1_en): (i32, i32, i32, i32)) {
        let d_en = (outer1_en + inner1_en) - (outer0_en + inner0_en);
        let a_en = (outer1_en - outer0_en).max(inner1_en - inner0_en);
        let bm = match (_o.is_stack(), _i.is_stack(), outer1.is_stack(), inner1.is_stack()) 
        {
            (true, false, false, true) => "neutral_trans",
            (false, true, true, false) => "neutral_trans",
            (false, true, false, true) => "neutral_cis", // loop
            (true, false, false, false) => "del_stack",
            (false, true, false, false) => "del_stack",
            (false, false, true, false) => "add_stack",
            (false, false, false, true) => "add_stack",
            (false, false, false, false) => "other",
            _ => panic!("strange rearrangement"),
        };
        println!("tw {} {} {} {}", d_en, a_en, 0.max(d_en).max(a_en), bm);
}

/// For a figure in Badelt et al. 2026: the energies of a four-way shift and
/// the kind of rearrangement, printed.
#[cfg(feature = "shift_analysis")]
fn print_four_way(_m1: &NestedLoop, _m2: &NestedLoop, outer1: &NestedLoop, outer2: &NestedLoop,
    (center_en, merge1_en, merge2_en, inside_en, outer1_en, outer2_en): (i32, i32, i32, i32, i32, i32)) {
        let d_en = (inside_en + outer1_en + outer2_en) - (center_en + merge1_en + merge2_en);
        let a_en = (outer1_en + outer2_en - center_en).max(inside_en - merge1_en - merge2_en);
        let bm = match (_m1.is_stack(), _m2.is_stack(), outer1.is_stack(), outer2.is_stack())
        {
            (true, true, true, true) => "neutral_2s",
            (true, false, true, false) => "neutral_1s",
            (true, false, false, true) => "neutral_1s",
            (false, true, true, false) => "neutral_1s",
            (false, true, false, true) => "neutral_1s",

            (true, false, true, true) => "add_stack",
            (false, true, true, true) => "add_stack",
            (false, false, true, false) => "add_stack",
            (false, false, false, true) => "add_stack",
            (false, false, true, true) => "add_stack",

            (true, false, false, false) => "del_stack",
            (false, true, false, false) => "del_stack",
            (true, true, true, false) => "del_stack",
            (true, true, false, true) => "del_stack",
            (true, true, false, false) => "del_stack",
            _ => "other",
            //(false, false, false, false) => "other",
            //_ => panic!("strange rearrangement"),
        };
        println!("fw {} {} {} {}", d_en, a_en, 0.max(d_en).max(a_en), bm);
}

#[cfg(test)]
mod tests {
    use super::*;
    use ff_energy::Fuzzyfold;
    use crate::movesets::loop_table::table;

    #[test]
    fn test_three_way_setup() {
        let ltab = table(Fuzzyfold::default(), "GUACACGUCG", "..........");
        let res = three_way_shifts(&ltab, 0);
        assert!(res.is_empty(), "Expected no three-way shift neighbors.");
    }
 
    #[test]
    fn test_three_way_simple() {
        let ltab = table(Fuzzyfold::default(), "GUACCCCUCG", "...(.....)");

        let res = three_way_shifts(&ltab, 0);
        assert!(res.len() == 2);
        let m1 = Move::ShiftJK { i: 3, j: 9, k: 4 };
        let m2 = Move::ShiftJK { i: 3, j: 9, k: 5 };
        assert_eq!(res[0].0, m1);
        assert_eq!(res[1].0, m2);

        let res = three_way_shifts(&ltab, 1);
        assert!(res.len() == 1);
        let m = Move::ShiftJK { i: 3, j: 9, k: 1 };
        assert_eq!(res[0].0, m);
    }
  
    #[test]
    fn test_three_way_mini() {
        let ltab = table(Fuzzyfold::default(), "CCCCUCG", "(.....)");

        let res = three_way_shifts(&ltab, 0);
        assert!(res.len() == 2);
        let m1 = Move::ShiftJK { i: 0, j: 6, k: 1 };
        let m2 = Move::ShiftJK { i: 0, j: 6, k: 2 };
        assert_eq!(res[0].0, m1);
        assert_eq!(res[1].0, m2);

        let res = three_way_shifts(&ltab, 1);
        assert!(res.is_empty());
    }

    #[test]
    fn test_four_way_setup() {
        let ltab = table(Fuzzyfold::default(), "GUACACGUCG", "..........");
        let res = four_way_shifts(&ltab, 0);
        assert!(res.is_empty(), "Expected no three-way shift neighbors.");
    }
 
    #[test]
    fn test_four_way_simple() {
        let ltab = table(Fuzzyfold::default(), "AGAAAAACAAAGAAACAA", ".(.....)...(...)..");

        let res = four_way_shifts(&ltab, 0);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");

        let res = four_way_shifts(&ltab, 1);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
        
        let res = four_way_shifts(&ltab, 2);
        let m = Move::ShiftILJK { i: 1, j: 7, k: 11, l: 15 };
        assert_eq!(res[0].0, m);

        let (_, _, _, inner, _, _) = four_way_loops(&ltab, &m);
        // The interior loop closed by (1, 15) with the inner pair (7, 11).
        assert_eq!((inner.gap(), inner.pairs()), (ff_energy::Gap::None, &[(7, 11), (15, 1)][..]));
    }

    #[test]
    fn test_four_way_simple_rev() {
        let ltab = table(Fuzzyfold::default(), "AGAAAAACAAAGAAACAA", ".(.....(...)...)..");

        let res = four_way_shifts(&ltab, 0);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
        
        let res = four_way_shifts(&ltab, 1);
        let m = Move::ShiftIKLJ { i: 1, k: 7, l: 11, j: 15 };
        assert_eq!(res[0].0, m);

        let res = four_way_shifts(&ltab, 2);
        assert!(res.is_empty(), "Expected no four-way shift neighbors.");
    }
}
