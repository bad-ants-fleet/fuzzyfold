# Milestones

## Train the pseudoknot parameters

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

Circular strands (evaluated by `Fuzzyfold`, also with pseudoknots and
other strands) raise two questions for the training:

- A circle joined to a linear strand is a pseudoknot: its ring spans the
  helix that joins them, so it pays the pseudoknot initiation, Pup per base
  of the ring and scaled stacks. Topologically right (a long duplex with a
  short circle is constrained), but with the untrained `PKParams` a short
  duplex on a long circle costs several kcal/mol more than on a linear
  strand. Compare with data.
- Reference energy of the open circle (no pairs): `Fuzzyfold` uses the
  hairpin initiation of its ring (n - 2 unpaired bases), `ViennaRNA` 0 (as
  ViennaRNA 2.6, which charges the ring only once a pair forms). Compare
  with experimental data on circular RNAs.

## Kinetics with pseudoknots and multiple strands

- Moves on pseudoknots, shift moves for all loop types. Relax the check in `LoopTable::try_from` as
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

## ff-findpath

A command line tool for the energy barrier between two structures A and B,
on `MoveSet` and the neighborhood enumeration (`enum_neighbors`), so it
extends to pseudoknots with a pseudoknot move set. Output: the saddle
energy and the path. Replaces the old `crates/ff_findpath` (removed; in
the git history).

- Direct paths (only pairs of A or B, every move one base pair closer to
  B): a breadth-first search by distance that keeps, for each structure,
  the lowest saddle of the paths to it.
  - `--fixed-width w`: keep the w structures with the lowest saddle per
    distance (findpath, Flamm et al. 2001; w = 1 is greedy). Without it,
    exhaustive: the exact direct barrier, exponential in the distance.
  - `--recursive`: take the estimated barrier as an upper bound for a new
    run, until the barrier does not change.
- `--indirect`: the full search over all paths of the move set, raising
  the energy bound until B is reached.
- `--stochastic`: a simple biased stochastic walk: with probability beta an
  SSA step, otherwise a step toward B (the biased walk of pathway
  elaboration, Zolaktaf et al. 2021).
- Open: shift moves on direct paths; restrict the move set to the allowed
  pairs instead of filtering the proposed moves; B as a macrostate.

## A fast heuristic for starting structures

A fast heuristic (pseudoknots allowed) for structures near dominant local
minima, expanded to macrostates with `ff-explore`, as starting points for
simulations of the equilibrium distribution.
