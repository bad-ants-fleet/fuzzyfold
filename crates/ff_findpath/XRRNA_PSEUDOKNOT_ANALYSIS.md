# xrRNA, Pseudoknot Energy Models, and Folding Paths — Scientific Notes

_Context: ff_findpath pseudo_findpath implementation, Rastegari-Condon loop enumeration,
Dirks-Pierce/Andronescu (DP09/MT09) pseudoknot energy model._

---

## 1. xrRNA and 7U4A — What We're Actually Looking At

### The xrRNA mechanism

xrRNAs (exoribonuclease-resistant RNAs) are structured elements in flavivirus genomes
(Zika, Dengue, West Nile, Japanese encephalitis) that stall the cytoplasmic 5'→3'
exonuclease XRN1 during degradation of viral RNA. The resistance is not energetic
(the structure can be thermally denatured) — it is **topological**: the RNA folds into
a ring-like geometry where the 5' terminus is threaded through or occluded by the
structure in such a way that XRN1 cannot processively unwind it.

The functional element is compact (~60–75 nt) and sits immediately downstream of a
cleavage site produced by a host endonuclease (e.g. Ago2 for some, SLIV cleavage for
others). After cleavage, the xrRNA fragment accumulates as a stable, non-coding sfRNA
(subgenomic flavivirus RNA).

### Structural topology of the xrRNA

The xrRNA typically contains:

1. **Two crossing helices** forming a pseudoknot at the 5' end of the element.
   - Helix P1 (5'-proximal): pairs with a region ≥10 nt downstream of the 5' terminus.
   - Helix P2 (3'-proximal): crosses P1 to pair with sequence between the P1 strands.
   - Together P1 + P2 constitute an H-type pseudoknot (two interleaved stems).

2. **Stem-loops** (P3, P4 in different nomenclatures) that dock into the major groove
   of P1 via tertiary contacts — notably, a conserved adenosine platform / base-triple
   interaction in the groove of P1. These are *not* Watson-Crick pairs and are invisible
   to secondary structure models.

3. **The ring**: The topology forces the 5' terminus to lie "inside" the ring defined
   by P1+P2 and the connecting loops. XRN1 cannot thread through this ring; it stalls
   at the 5' face of P1.

### On PDB 7U4A

7U4A is a cryo-EM or X-ray structure of a flaviviral xrRNA deposited in the RCSB PDB.
The specific organism / publication tied to 7U4A should be verified at rcsb.org before
drawing structure-specific conclusions — this analysis proceeds at the level of xrRNA
class features, which are well-conserved.

Key structural data to look at in 7U4A:
- Which positions form P1 vs. P2 (extract pair table from ATOM records)
- Whether the pseudoknot is a classical H-type or involves a third crossing helix
  (some xrRNAs have three bands, not two — the DP09 model only handles two reliably)
- The loop lengths between helices: these set the `n_loop1`, `n_loop2` parameters
  in our pseudoloop energy term
- Whether any non-Watson-Crick base triples are present in the loops
  (these contribute to stability but are absent from our energy model)

**Practical step**: extract the secondary structure from 7U4A ATOM records with
`x3dna-dssr` or RNAView, convert to dot-bracket, run through `parse_structure` to get
the loop decomposition and pseudoloop energy term we'd assign to this structure.
Compare to experimental Tm data if available.

---

## 2. What Does the findpath Path Tell Us About the PK Energy Model?

### What findpath does (and assumes)

`findpath_pseudo` finds the minimum-saddle-energy directed path between two structures.
"Directed" means only moves that make progress toward the target are considered — there
is no backtracking, no alternative misfolded intermediate exploration.
The energy at each intermediate is evaluated using the DP09/MT09 parameter set via
the Rastegari-Condon loop enumeration.

### What the path actually encodes about the model

**2a. The per-band penalty, not the initiation penalty, drives the visible energy spike.**

The full DP09 pseudoloop junction formula is:

```
G_pseudoloop = init(ctx) + pb × n_bands + pup × (n_loop1 + n_loop2) + pps × n_nested
             + e_stp  × Σ stacking_energies_in_bands
             + e_intp × Σ interior_loop_energies_in_bands
```

A common simplification in the literature writes this as `a_pk + b_pk × n_unpaired`,
where `a_pk` is described as "a large positive constant (~9 kcal/mol)". That figure
comes from **dp03** (`init_external = +960 dcal/mol = +9.60 kcal/mol`). In the
**dp09** parameter set used here, the cost is redistributed:

| Parameter | dp03 | dp09 |
|---|---|---|
| `init_external` (per pseudoloop, exterior context) | +9.60 kcal/mol | **−1.38 kcal/mol** |
| `pb` (per crossing helix / band) | +0.20 kcal/mol | **+2.46 kcal/mol** |
| `e_stp` (stack energy scale inside band) | 0.83 | 0.89 |

For `[[[....(((((]]]......)))))` the pseudoknot is in exterior context and has
**2 bands**. The moment the first crossing pair is inserted, both bands are present
and the full per-band penalty is paid:

```
pb × n_bands = 2.46 × 2 = +4.92 kcal/mol   (dp09)
```

Simultaneously, all stacking energies inside the now-crossing helices are scaled by
`e_stp = 0.89` — an 11% reduction applied to whatever helix was already there.

This explains the sharp energy spike visible at the nucleation step even though
`init_external` is favourable (−1.38 kcal/mol): in dp09 the "cost" of forming a
pseudoknot was not eliminated, it was **redistributed** from a single large initiation
penalty into per-band and scaling penalties. The two models predict similar nucleation
barriers by design (both were fitted to the same experimental data); only the
parameter interpretation differs.

Consequence for the findpath path:
- There will be a **sharp energy spike** at the step where the first crossing pair
  is inserted, because that step pays `pb × n_bands` in full plus the `e_stp` scaling
  loss on any pre-existing helix — without yet gaining stacking from the new helix.
- For H2-preformed → PK, the spike is especially large because H2 is very stable
  (−6.92 kcal/mol) and the `e_stp` scaling loss on five already-formed stacks is
  significant (+0.7 kcal/mol at 0.89 scale).
- The path will tell us whether forming **P1 first vs P2 first** gives a lower barrier —
  this is a prediction of the energy model about the preferred nucleation order.

**2b. The pseudoloop junction parameters are the most uncertain.**

The DP09 pseudoloop parameters (`init_external`, `pb`, `pup`, `pps`, `e_stp`,
`e_intp`) were not fit to folding-kinetic data but to equilibrium UV-melting curves
of pseudoknotted RNAs. The parameter covariance is high: `init_external` and `pb`
are strongly anti-correlated (lowering initiation and raising per-band penalty can
give equally good fits to the training data). Different fitting runs give total
junction penalties ranging from roughly +3 to +12 kcal/mol, even when each
individual fit is thermodynamically valid.

This means that the absolute barrier heights from findpath for PK-forming steps have
large systematic uncertainty — easily ±3 kcal/mol.

**What this means for interpretation:**
- The *relative* ordering of paths (e.g. P1-first vs P2-first barrier comparison) is
  more trustworthy than absolute values, because the full junction penalty appears in
  both and partially cancels.
- The stacking contributions inside each helix are well-parameterised (Turner 2004
  basis); the uncertainty is concentrated in the junction parameters.

**2c. The path reveals energy landscape topology, not just the minimum path.**

At each step, the saddle energy recorded by findpath is the minimum achievable saddle
*over all orderings of the remaining moves*. Inspecting the saddle progression reveals:

- If the saddle is reached early (first few moves), the model sees a high barrier to
  nucleation that is not subsequently exceeded. This is typical for H-type PK formation:
  once you pay `a_pk`, subsequent insertions add stacking that is stabilizing.
- If the saddle is reached late (last few moves), the model says one of the final
  insertions is the highest-energy step. This would be unusual for a cooperative helix
  and might indicate a steric clash or mispairing cost in the model.
- A flat saddle profile (all energies near the saddle) suggests the model sees the
  transition state as broad — the barrier is not a single step but a plateau.

**Concrete test to run:**
Compare the barrier for (empty → xrRNA) against:
1. (empty → P1 only) + (P1-only → xrRNA) — does P1 form as a stepping stone?
2. (empty → P2 only) + (P2-only → xrRNA) — does P2 form as a stepping stone?

If the barrier for path 2 is significantly lower, the model predicts P2 nucleates first
despite the 5'-distal position of P2 in the sequence. This would be a testable
prediction against SHAPE-MaP or NMR time-resolved experiments.

---

## 3. Is the PK Energy Model Useful for Cotranscriptional Folding Paths?

### The cotranscriptional folding problem

During transcription, the RNA is synthesized 5'→3'. At each elongation step, one new
nucleotide is added and can pair with any upstream position. The RNA does not wait
until synthesis is complete before folding — each new segment can immediately form
secondary (and potentially tertiary) structure.

For xrRNA, cotranscriptional folding is biologically critical: the structure must
assemble as the virus RNA is synthesized (or as the host ribosome reads it), before
the xrRNA region is completely transcribed.

### What findpath can contribute to this picture

findpath is not cotranscriptional by design, but it can be adapted:

**Approach A: Transcript-length series**

For each transcript length L from `L_min` (where the first PK pair becomes possible)
to `L_full` (complete xrRNA), run findpath from the empty structure at length L to
the optimal structure of the length-L prefix. The barriers at each length reveal when
the PK-forming step first becomes accessible and whether it is kinetically favoured
at that length.

This is well within what `findpath_pseudo` supports via the `start: Option<&str>`
parameter — at each length, the start can be the winner from the previous length.

**Approach B: Non-empty start as proxy for cotranscriptional state**

The most direct use of the non-empty start feature: set the start to the structure
adopted by the partially synthesized transcript (from a ViennaRNA MFE calculation on
the L-nt prefix), and find the path to the full xrRNA. The barrier in this direction
is the kinetic cost of "switching" from the naive co-transcriptional fold to the
functional pseudoknot.

If this barrier is high (> ~5 kcal/mol at 37°C), the model predicts kinetic trapping
in non-pseudoknotted intermediates. If it is low, the model predicts smooth
cotranscriptional assembly.

### Fundamental limitations of the DP09 model for this purpose

**L1. The junction penalty does not depend on transcript length.**

The `init(ctx) + pb × n_bands + …` terms are constants given the loop topology. In
reality, the barrier to form a pseudoknot on a growing chain depends on whether both
strands of P1 and P2 are already synthesized. Before the 3′-strand of P2 is
transcribed, P2 cannot form — but the model would still assign the per-band penalty
for a structure with only P1 present. The model handles missing pairs implicitly
(unpaired positions), but the *pathway* through partial pseudoknots during synthesis
is not captured by a static energy function.

**L2. The model only handles H-type pseudoknots (2 bands) robustly.**

More complex topologies that appear transiently during cotranscriptional folding — for
instance, a structure where P1 is formed but P2 has only 1 pair and also overlaps with
a stem-loop from the upstream sequence — may not be H-type and will either fail energy
evaluation or be assigned incorrect energy. This is the `Err(_) => continue` we see in
the beam search: the model silently skips topologically complex intermediates.

This is a real scientific blind spot: the intermediates that are most interesting for
cotranscriptional folding may be precisely those that are topologically complex.

**L3. Tertiary contacts are invisible.**

The base-triple interactions and loop-groove docking that stabilize the xrRNA structure
contribute ~2–5 kcal/mol to stability (estimated from mutant thermodynamics) but are
entirely absent from the nearest-neighbor model. The findpath energy landscape is
therefore smoother and less funnelled than the actual landscape.

**L4. The model was not trained on kinetic data.**

The parameters are thermodynamic (equilibrium) values. The Arrhenius rates used in
ff_kinetics (`k0`, `k3ws`, `k4ws`) are physically motivated but not calibrated to
xrRNA-specific kinetics. Absolute rate predictions should not be taken at face value.

### Where the model is nonetheless useful

Despite limitations L1–L4, the DP09+MT09 model is the only available analytical model
that:
- Assigns a meaningful energy to pseudoknotted RNA structures (unlike the ViennaRNA
  2-structure model which cannot evaluate PKs at all)
- Has a physically grounded loop-decomposition scheme (Rastegari-Condon) that maps
  to identifiable structural elements
- Is fast enough to evaluate thousands of intermediates (especially after the B1
  optimization: ~17 ms for 70 nt / beam=10)

The scientifically valid uses are:

1. **Comparative analysis**: Compare the barrier to reach xrRNA structure from empty
   vs. from a misfolded intermediate. The *ratio* of barriers is less sensitive to
   parameter uncertainty than absolute values.

2. **Nucleation order prediction**: The model robustly identifies which helix (P1 or
   P2) is the preferred nucleation site, because this is dominated by stacking
   energetics (well-parameterised) relative to each other. The junction penalty
   (`pb × n_bands`, `init_external`) is the same for both nucleation orders and
   largely cancels in the comparison.

3. **Identifying kinetic traps**: If a misfolded intermediate has a high barrier to
   escape (> 3–4 kcal/mol), the model is likely correct that it represents a kinetic
   trap, even if the exact value is off. This is because the barrier is the *difference*
   in energy between a stable stem and the partially disrupted transition state —
   the stacking terms dominate, and these are well-parameterised.

4. **Screening**: Run findpath across many xrRNA variants (natural or mutant) to
   screen for variants with unusually high or low PK formation barriers. Relative
   rankings are more reliable than absolute barriers.

---

## Summary Table

| Question | What we can say | Uncertainty |
|----------|----------------|-------------|
| Does xrRNA fold through a defined intermediate? | Yes — model predicts P1- or P2-nucleation as a preferred intermediate; findpath identifies which | ±3 kcal/mol on barrier height from junction parameter uncertainty (`init_external`, `pb` covariance) |
| Which helix forms first? | findpath gives a testable prediction (P1-first vs P2-first path barrier comparison) | Stacking terms reliable; initiation penalty partially cancels in comparison |
| Is the PK model useful for cotranscriptional paths? | Yes for relative comparisons and kinetic trap identification; no for absolute rate or barrier prediction | Tertiary contacts (2–5 kcal/mol) and complex topologies are systematically missing |
| Can we use non-empty start for cotranscriptional simulation? | Yes — transcript-length series with start = MFE of prefix is a principled approach | Intermediate structures may be topologically complex and evaluator will silently skip them |

---

## Suggested Next Steps

1. **Extract pair table from 7U4A** using DSSR or similar; run `parse_structure` to
   get the Rastegari-Condon loop decomposition; note the `n_loop1`, `n_loop2` values
   and the number of bands (if > 2, the model needs extension).

2. **Run findpath_pseudo on the 7U4A sequence**:
   - empty → xrRNA structure (barrier A)
   - MFE of full sequence without PK → xrRNA structure (barrier B)
   - If B >> A, the model predicts misfolding is a problem.

3. **Transcript-length series**: at each L from first-PK-possible to full length,
   compute the findpath barrier for the PK-forming move. Plot barrier vs. L to see
   the window where PK formation becomes kinetically accessible.

4. **Compare P1-first vs P2-first nucleation** by running:
   - findpath(empty → P1-only) then findpath(P1-only → full PK)  
   - findpath(empty → P2-only) then findpath(P2-only → full PK)
   - Sum each and compare total barriers.

5. **Quantify tertiary-contact gap**: If you have Tm data for the xrRNA, compare the
   experimental stability to the DP09-predicted `ΔG_37`. The gap is the tertiary
   contribution missing from the model — this gives an upper bound on how much the
   barrier could be underestimated.
