//! [`PseudoEnergyModel`] implementation for [`ViennaRNA`].

use ff_structure::NAIDX;

use super::ViennaRNA;
use crate::pseudoknots::{
    ClosingDescriptor, LocationStatus, Loop, LoopType, PseudoEnergyModel, PseudoloopContext,
    double_pair, single_pair,
};
use crate::{Base, EnergyError, EnergyModel, NearestNeighborLoop};

impl PseudoEnergyModel for ViennaRNA {
    fn energy_of_pseudo_loop(&self, sequence: &[Base], lp: &Loop) -> Result<i32, EnergyError> {
        match lp.loop_type {
            // ── Step 6: SpanBand stacks / interior loops use scaled Turner energies ──
            LoopType::Stack | LoopType::Interior | LoopType::Bulge => {
                let (i, j) = single_pair(lp.closing)?;
                let (k, l) = lp.inner.ok_or(EnergyError::InvalidClosingPair)?;
                let raw = self.energy_of_loop(
                    sequence,
                    &NearestNeighborLoop::Interior {
                        closing: (i as NAIDX, j as NAIDX),
                        inner: (k as NAIDX, l as NAIDX),
                    },
                )?;
                if lp.location == LocationStatus::SpanBand
                    && let Some(pk) = &self.pk_params
                {
                    let scale = if lp.loop_type == LoopType::Stack {
                        pk.e_stp
                    } else {
                        pk.e_intp
                    };
                    return Ok((raw as f64 * scale).round() as i32);
                }
                Ok(raw)
            }

            LoopType::Hairpin => {
                let (i, j) = single_pair(lp.closing)?;
                self.energy_of_loop(
                    sequence,
                    &NearestNeighborLoop::Hairpin {
                        closing: (i as NAIDX, j as NAIDX),
                    },
                )
            }

            // ── Step 5: SpanBand multiloops use ap/bp/cp instead of Turner ──
            LoopType::Multiloop => {
                if lp.location == LocationStatus::SpanBand
                    && let Some(pk) = &self.pk_params
                {
                    let n_branches = lp.children.len();
                    let n_unpaired = lp.unpaired_5p + lp.unpaired_3p;
                    return Ok(pk.ap + pk.bp * n_branches as i32 + pk.cp * n_unpaired as i32);
                }

                let (i, j) = single_pair(lp.closing)?;
                // A pseudoknot child has a Double closing descriptor with crossing
                // pairs (a1 < a2 < b1 < b2). We use the leftmost pair (a1, b1)
                // as the effective multiloop stem: it is a real WC pair, so the
                // pair checks in `multibranch` pass. The region b1..b2 (the far
                // strand of the PK) is counted as loop nucleotides — a small
                // approximation. The PK's internal energy comes from its own
                // LoopType::Pseudoloop entry.
                let branches: Vec<(NAIDX, NAIDX)> = lp
                    .children
                    .iter()
                    .map(|cd| match cd {
                        ClosingDescriptor::Single(p) => (p.0 as NAIDX, p.1 as NAIDX),
                        ClosingDescriptor::Double(p1, _p2) => (p1.0 as NAIDX, p1.1 as NAIDX),
                    })
                    .collect();
                self.energy_of_loop(
                    sequence,
                    &NearestNeighborLoop::Multibranch {
                        closing: (i as NAIDX, j as NAIDX),
                        branches,
                    },
                )
            }

            LoopType::External => {
                // Each branch (Single or Double descriptor) contributes independently:
                // AU/GU terminal penalty + dangling ends based on adjacent sequence context.
                // Crossing pairs from a Double descriptor cannot share a single slice pass,
                // so every pair is evaluated by slicing [..=i] and [j..] individually.
                let mut en = 0;
                for cd in &lp.children {
                    match cd {
                        ClosingDescriptor::Single((a, b)) => {
                            en += self.exterior(&[&sequence[..=*a], &sequence[*b..]])?;
                        }
                        ClosingDescriptor::Double((a1, b1), (a2, b2)) => {
                            en += self.exterior(&[&sequence[..=*a1], &sequence[*b1..]])?;
                            en += self.exterior(&[&sequence[..=*a2], &sequence[*b2..]])?;
                        }
                    }
                }
                Ok(en)
            }

            // ── Steps 1-4: Full 6-term D&P Pseudoloop formula ──
            LoopType::Pseudoloop => {
                let _ = double_pair(lp.closing)?; // validate descriptor
                if let Some(pk) = &self.pk_params {
                    // E = init(ctx) + pb×n_bands + pup×n_pup + pps×n_nested
                    let init = match lp.pk_context.unwrap_or(PseudoloopContext::External) {
                        PseudoloopContext::External => pk.init_external,
                        PseudoloopContext::Multiloop => pk.init_multiloop,
                        PseudoloopContext::Pseudoloop => pk.init_pseudoloop,
                    };
                    Ok(init
                        + pk.pb * lp.n_bands as i32
                        + pk.pup * lp.n_pup as i32
                        + pk.pps * lp.n_nested as i32)
                } else {
                    // Fallback: simplified multiloop approximation.
                    let n_stems = lp.children.len() + 2;
                    Ok(self.ml_closing + self.ml_intern * n_stems as i32)
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameters::RNA_EXTENDED;
    use crate::{NucleotideVec, parse_structure};

    fn model() -> ViennaRNA {
        ViennaRNA::from_thermo_params(&RNA_EXTENDED, 37.0)
    }

    /// A non-pseudoknotted structure should give the same total energy via
    /// `energy_of_pseudoknotted_structure` and the standard `energy_of_structure`.
    #[test]
    fn test_non_pk_matches_standard() {
        use ff_structure::PairTable;

        let seq = NucleotideVec::try_from_rna("GGGAAACCC").unwrap();
        let dot = "(((...)))";
        let pt = PairTable::try_from(dot).unwrap();
        let loops = parse_structure(dot).unwrap();

        let m = model();
        let e_standard = m.energy_of_structure(&seq, &pt).unwrap();
        let e_pseudo = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert_eq!(
            e_standard, e_pseudo,
            "non-PK energy mismatch: standard={e_standard}, pseudo={e_pseudo}"
        );
    }

    /// An H-type pseudoknot should evaluate without error.
    #[test]
    fn test_h_type_pseudoknot_runs() {
        // GGGCCCCCCAAAGGG pairs as ((([[[)))...]]]
        let seq = NucleotideVec::try_from_rna("GGGCCCCCCAAAGGG").unwrap();
        let loops = parse_structure("((([[[)))...]]]").unwrap();

        let m = model();
        let result = m.energy_of_pseudoknotted_structure(&seq, &loops);
        println!("Pseudoknot energy: {:?}", result);
        assert!(
            result.is_ok(),
            "pseudoknot energy evaluation failed: {:?}",
            result
        );
    }

    /// Smoke test: dp09 junction formula is pure arithmetic, no Turner tables needed.
    ///
    /// Same structure as above.
    ///   init_ext=−138 + pb=246×2 + pup=6×(4+6) + pps=96×0 = −138+492+60 = 414
    #[test]
    fn test_dp09_junction_formula() {
        use crate::parameters::{RNA_DP09, RNA_TURNER_2004};

        let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("[[[....(((((]]]......)))))").unwrap();
        let m =
            ViennaRNA::from_thermo_params(&RNA_TURNER_2004, 37.0).with_pseudoknot_params(RNA_DP09);
        let pk_loop = loops
            .iter()
            .find(|l| l.loop_type == LoopType::Pseudoloop)
            .unwrap();

        assert_eq!(m.energy_of_pseudo_loop(&seq, pk_loop).unwrap(), 414);
    }

    /// Full dp09+mt09 total energy for an H-type pseudoknot.
    ///
    /// HotKnots v2 reference (with mt09 NN + dp09 PK): −5.11 kcal/mol = −511 dcal/mol.
    /// dp09 PK parameters were trained jointly with the mt09 NN parameters; using them
    /// with Turner 2004 tables instead gives systematically wrong results (~−3.43 kcal/mol).
    /// Tolerance ±1 dcal/mol for f64→i32 rounding in stack scaling.
    #[test]
    fn test_dp09_total_energy() {
        use crate::parameters::{RNA_DP09, RNA_MT09};

        let seq = NucleotideVec::try_from_rna("GCGAUUUCUGACCGCUUUUUUGUCAG").unwrap();
        let loops = parse_structure("[[[....(((((]]]......)))))").unwrap();
        let m = ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09);
        let e = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert!(
            (e - (-511)).abs() <= 1,
            "dp09+mt09 total energy: expected ~−511 dcal/mol, got {e}"
        );
    }

    /// A four-band structure containing a nested pseudoknot: exercises the
    /// `n>2` band path and the `init_pseudoloop` context together.
    ///
    /// No external reference exists for this one — dp09 was trained on 2-band
    /// structures, so a 4-band value is an extrapolation. The assertion is an
    /// internal regression baseline: it pins the value so that a change in the
    /// enumeration or the DP terms shows up rather than passing silently.
    #[test]
    fn test_dp09_four_band_nested_pk() {
        use crate::parameters::{RNA_DP09, RNA_MT09};

        let seq = NucleotideVec::try_from_rna(
            "ACCCCAAACCCCAAAGGGGAAACCCCAAAGGGGAAACCCCAAAGGGGAAACCCCAAACCCCAAAGGGGAAAGGGGAGGGGA",
        )
        .unwrap();
        let loops = parse_structure(
            ".((((...[[[[...))))...((((...]]]]...[[[[...))))...((((...[[[[...))))...]]]].]]]].",
        )
        .unwrap();

        // Structural preconditions the energy value depends on.
        let outer = loops
            .iter()
            .find(|l| l.loop_type == LoopType::Pseudoloop && l.n_bands == 4)
            .expect("outer 4-band pseudoloop");
        assert_eq!(outer.pk_context, Some(PseudoloopContext::External));
        assert_eq!(outer.n_nested, 1, "contains one nested pseudoknot");

        let inner = loops
            .iter()
            .find(|l| l.pk_context == Some(PseudoloopContext::Pseudoloop))
            .expect("nested pseudoloop");
        assert_eq!(inner.n_bands, 2);

        let m = ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09);
        let e = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        assert_eq!(e, -446, "internal regression baseline (dp09+mt09)");
    }

    /// Zika xrRNA1 — a real biological pseudoknot, four bands with a nested
    /// pseudoknot inside. Internal regression baseline, as above.
    #[test]
    fn test_dp09_zika_xr_rna1() {
        use crate::parameters::{RNA_DP09, RNA_MT09};

        let seq = NucleotideVec::try_from_rna(
            "AAGUCGGCCGUACGAAUAGUACCCAGGACAGCCUCCCUCUGUCCUGCUGGCCGUGGGAGAG",
        )
        .unwrap();
        let loops =
            parse_structure(".[[.(((((((((.....)))).((((((((.<<<<<.))))))))]]))))).>>>>>..")
                .unwrap();

        let pk = loops
            .iter()
            .find(|l| l.loop_type == LoopType::Pseudoloop)
            .expect("pseudoloop");
        assert_eq!(pk.n_bands, 4);

        let m = ViennaRNA::from_andrunescu_params(&RNA_MT09).with_pseudoknot_params(RNA_DP09);
        let e = m.energy_of_pseudoknotted_structure(&seq, &loops).unwrap();
        // Baseline moved -1586 -> -1510 -> -1693 across three fm363 converter
        // fixes: the hairpin-block loop bound, folding the tetraloop bonus into
        // an absolute energy, and building MISMATCH_MULTI / MISMATCH_EXTERIOR
        // from the dangle model instead of copying MISMATCH_INTERIOR. The
        // external HotKnots reference in test_dp09_total_energy never moved.
        assert_eq!(e, -1693, "internal regression baseline (dp09+mt09)");
    }
}
