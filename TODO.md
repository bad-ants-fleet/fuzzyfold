# Milestones

## 1. Merge pseudoknots and multiple strands (energy evaluation)

Finish the review of this branch and merge it into `development`, then `main`.

- Review `ff_energy`: `models/fuzzyfold.rs` (`Fuzzyfold::pseudoknot()`),
  `decomposition/pseudoknots.rs` (the loop graph, strained segments,
  knots) and `decomposition/rings.rs`.
- Against main (instruction counts, same parameters, identical results):
  evaluation about 10% faster; kinetics about 12% (no shifts) and 17%
  (shifts) slower, mostly from loop energies split into separate functions
  (`dangles`) and segments built for multiloops. Measure before changing.
- Crate docs for `ff_kinetics`, like the guides in `ff_structure` and
  `ff_energy`; examples for docs.rs.
- CHANGELOGs of `ff_kinetics` and the Python bindings (Python: `ff.Fuzzyfold`
  and `ff.ViennaRNA`, default parameters `rna_extended`); version bump.
- Citation: `CITATION.cff`, a citation line in the `--help` of the command
  line tools, `pyproject.toml` metadata (readme, authors, URLs, license).
- Decide the fate of `crates/ff_findpath` (not in the workspace, outdated).

## 2. Train the pseudoknot parameters

The parameters belong to `Fuzzyfold`, the model under development.
`ViennaRNA` stays the model of ViennaRNA 2.6, without pseudoknots, for
comparisons. `Fuzzyfold` uses the loop energies of `ViennaRNA` until it
diverges (e.g. a new multiloop model or coaxial stacking): then it gets its
own function for that loop type, also where `Fuzzyfold::pseudoknot()`
evaluates regular loops inside a pseudoknot.

Before any release with pseudoknots: the defaults come from a preliminary
training on the HotKnots data sets.

- An Andronescu-style interface: target sequence-structure pairs, some with
  known energies; vary only selected parameters (pseudoknot, multiloop).
- Constraints: the target is not above any structure in its local
  neighborhood (hard); external knots (soft). Decoys come from the local
  neighborhood.
- Decide the loop model (current a + b * branches + c * unpaired, or a
  log-length model) before training, not because of it.
- One trained default per nearest neighbor set (Turner 2004, Andronescu
  2007). Temperature dependence of the pseudoknot terms (today: a warning).

## 3. Manuscript: multistranded pseudoknot evaluation

The method (stems, rings, strand order independence), comparison with
existing models (HotKnots / Dirks & Pierce), arc plot figures. State the
difference to HotKnots: Pb per stem, not per band (a band of Rastegari &
Condon may span a multiloop, a stem ends there), and no separate parameters
for a multiloop inside a band (a', b', c'). Parameter
training as a separate part that does not drive the model choices.

## 4. Next release: kinetics with pseudoknots and multiple strands

- Moves on pseudoknots, shift moves for all loop types (`todo!` in
  `three_way_shifts.rs`). Relax the check in `LoopTable::try_from` as
  support grows. `split_loop` / `join_loop` work on loops that meet their
  pairs in order only (`NestedLoop` of a nested structure). A pair added
  in a pseudoknot lies in an interior loop of a stem (splits
  the stem), between stem regions (a new child), encloses branches (they
  move inside), or crosses stem pairs (changes the stems); the last case,
  e.g. extending a stem, does not fit an (outer, inner) result. Removing a
  pair next to a pseudoknot may shrink or remove a stem, or dissolve the
  pseudoknot into nested regions.
- Multiple strands: joint exterior loops, association and dissociation
  (rates depend on concentration), detached complexes.
- Macrostates with pseudoknots: a `DotBracketVec` key needs canonical bracket
  levels at every step; compare with `MultiPairTable` keys.
- Runtime expectations: little existing work to compare against (Multistrand
  for multiple strands without pseudoknots). Set up our own benchmarks early.

## Circular strands

Represented (`^` in sequence and structure). Circular strands are
evaluated, also with pseudoknots and with other strands (`Fuzzyfold`;
`ViennaRNA` only a circle alone): the loop through the wrap is closed by
the backbone, so it is strained once it lies on a cycle, like a hairpin
loop that forms a kiss; a stem may continue through the wrap. Energies do
not depend on where the circle is cut (tests in `pseudoknot_energies.rs`
and `test_circular_rotation_invariance`). Simulated at its full length,
without pseudoknots (no co-transcriptional folding).

- A circle joined to a linear strand is a pseudoknot: its ring spans the
  helix that joins them, so it pays the pseudoknot initiation, Pup per base
  of the ring and scaled stacks. Topologically right (a long duplex with a
  short circle is constrained), but with the untrained `PKParams` a short
  duplex on a long circle costs several kcal/mol more than on a linear
  strand. Compare with data before training.
- Reference energy of the open circle (no pairs): `Fuzzyfold` uses the
  hairpin initiation of its ring (n - 2 unpaired bases), `ViennaRNA` 0 (as
  ViennaRNA 2.6, which charges the ring only once a pair forms). Compare
  with experimental data on circular RNAs before the paper.
- The pseudoknot density (`ClosedRegions`) depends on where the circle is
  cut.

## Later

- A `model=` argument for the Python `Simulator` and `Explorer`, like
  `--model` of the command line tools. The CLI enum `energy_parsers::Model`
  would move to where both can use it.
- `energy_terms` as a method of `EnergyModel` (today one per model).
- `andronescu2007` in the Python bindings.

- A fast heuristic (pseudoknots allowed) for structures near dominant local
  minima, expanded to macrostates with `ff-explore`, as starting points for
  simulations of the equilibrium distribution.
