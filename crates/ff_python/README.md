# fuzzyfold — Python bindings

Python bindings for fuzzyfold, a Rust framework for nucleic acid energy 
evaluation and stochastic folding kinetics. The bindings are a thin layer 
over the Rust crates: models, simulations and file formats are the same as 
in the `ff-*` command line tools.

## Installation

Requires a Rust toolchain and Python ≥ 3.10. From this directory:

```bash
pip install .                 # build and install the module
maturin develop --release     # or, during development, into the active environment
```

The version follows the crate (`Cargo.toml`).

All energies are in kcal/mol.

## Energy evaluation

```python
import fuzzyfold as ff

model = ff.Fuzzyfold()                        # rna_extended at 37 C
model = ff.Fuzzyfold(params="dna", celsius=25)
model = ff.ViennaRNA(celsius=25)              # as ViennaRNA 2.6, no pseudoknots
model.energy_of_structure("UGCCUAGAGAGUCAGGUGAU", ".((((.((...))))))...")   # kcal/mol
```

`params` is `"rna_extended"` (Turner 2004 with special hairpins and
modified bases; default of `Fuzzyfold`), `"rna_turner_2004"` (default of
`ViennaRNA`) or `"dna"` (Mathews 2004). 
Strands are separated by `+` or `&`. A `^` after a strand makes it circular,
in sequence and structure. `Fuzzyfold` evaluates circular strands also with
pseudoknots and together with other strands (`ViennaRNA` only a circle
alone); kinetics simulate a circular strand alone, without pseudoknots.
`ViennaRNA` warns for circular strands: it scores the open circle as 0,
which is not the reference state of the structured circles.

### Pseudoknots

```python
model = ff.Fuzzyfold()
seq, pk = "GCGAUUUCUGACCGCUUUUUUGUCAG", "[[[....(((((]]]......)))))"
model.energy_of_structure(seq, pk)            # -3.57
for level, term, energy in model.energy_terms(seq, pk):
    print("  " * level, term, energy)
ff.stems(pk)                                  # [(pairs, strained), ...]
```

Pseudoknots use `()[]{}<>` and the letter families `Aa` to `Zz`. They are 
evaluated with pseudoknot terms on top of the nearest neighbor model, the
same as `ff-eval`.
For now, every parameter set uses the same pseudoknot parameters, from a
preliminary training on the HotKnots data sets, not rescaled with
temperature.

- `energy_terms` lists each loop and, for a pseudoknot, its terms (the
  loops of its strained stems, and the helix ends of its loops), with their
  level: the content of `ff-eval -v`.
- `stems` returns the stems of the pseudoknots as the evaluation 
  sees them, each as its pairs and whether it is strained (scaled, Pb) or
  relaxed.

## Structures

Functions on dot-bracket strings. Indices are positions in the string, so a
strand break (`+`) or strand wrap (`^`) takes one position. Invalid
structures raise `ValueError`.

```python
s = "((..[[..))..]]"
ff.partners(s)          # [9, 8, None, None, 13, 12, None, None, 1, 0, ...]
ff.pairs(s)             # [(0, 9), (1, 8), (4, 13), (5, 12)]
ff.normalize("[(.)]")   # '((.))': one spelling per structure
ff.closed_regions(s)    # [(0, 13, False)]: (i, j, nested)
ff.density(s)           # 2
```

- `normalize` assigns the brackets again from `(` upward and writes strand
  breaks as `+`, so equal structures give equal strings, e.g. as dict keys.
- `closed_regions`: intervals that no pair crosses, minimal from their first
  base; a region is nested if its first and last base pair with each other.
- `density`: the largest number of bands of one pseudoknotted region that
  span a common position; 1 for nested structures, 2 for H-types and
  kissing hairpins.

## Kinetics

- `Simulator(params, celsius, k0, k3ws, k4ws, k_atoi)`: stochastic
  simulation (Gillespie), optionally co-transcriptional (`t_ext`), with
  three- and four-way branch migration (`k3ws`, `k4ws`), and with A-to-I
  editing (`k_atoi`, see below).
  - `simulate_trajectory(...)`: an iterator over the trajectory.
  - `simulate_ensemble(...)`: the structures of many simulations at the
    output times, also written as a `.drf` file.
  - `simulate_timecourse(...)`: the macrostate counts of many simulations,
    as a `Timecourse`; with `output`, also written as `.tln` and `.nxy`
    files (the files of `ff-timecourse`).
- `Explorer(...).explore(sequence, structure, ...)`: the neighborhood of a 
  structure within an energy range.

Pseudoknots are not supported in kinetics yet.

With `k_atoi > 0`, an A opposite a C in a 1x1 interior loop becomes
inosine at that rate, and the I-C pair forms. The sequence then changes
during a simulation: the trajectory iterator has the sequence of its last
state as `iterator.sequence`; `simulate_ensemble` counts each structure
together with its sequence (keys `"structure sequence"`, and a sequence
column in the `.drf` file); a macrostate file may give an edited sequence
(I for A), and a structure then belongs to the macrostate of its
sequence.

`simulate_ensemble` and `simulate_timecourse` report the state at the
output times `times`, and a simulation ends at the last one. By default,
they are 50 times evenly up to 1e-4 s and 50 in log space up to 1 s
(`ff.lin_log_times()`); with `t_ext`, one per extension up to the full
length and 50 in log space up to 1 s after it
(`ff.transcription_times(t_ext, num_ext)`). Any increasing list works:

```python
times = ff.lin_log_times(t_sep=1e-5, t_end=10, t_lin=20, t_log=100)
```

### Time courses

A `Timecourse` holds the counts of each macrostate (per sequence length
in co-transcriptional simulations) at the output times, and the number of
simulations. It needs numpy.

```python
sim = ff.Simulator(k0=1e5)
tc = sim.simulate_timecourse(seq, start, num_sims=100,
                             macrostates=["lm1.ms", "lm2.ms"])
tc.names, tc.times, tc.occupancy()      # occupancy: (times, names)
tc.standard_error()
tc.save("run.tln")                      # or .nxy (occupancies only)
tc = ff.Timecourse.read("run1.tln").merge(ff.Timecourse.read("run2.tln"))
```

### Plots

`fuzzyfold.plotting` plots time courses with matplotlib
(`pip install .[plot]`): linear up to `t_sep` and logarithmic after it,
with the standard errors as bands, and a sequence length axis for
co-transcriptional simulations. A `Timecourse` or a `.tln` file gives
`t_sep` and the errors; for a `.nxy` file, give `t_sep` and `num_sims`.

```python
from fuzzyfold.plotting import plot_timecourse
fig = plot_timecourse(tc, "run.pdf", title="100 simulations")
```

```bash
python -m fuzzyfold.plotting run1.tln run2.tln -o run -f pdf png
python -m fuzzyfold.plotting run.nxy --t-sep 1e-5 --labels LM1 LM2 -f show
```

`python -m fuzzyfold.plotting -h` lists the options.

## Examples

- `examples/example.py`: energy evaluation, simulations, and the energy 
  breakdown of a pseudoknot.
- `examples/ff_eval_arcs.py`: `ff-eval` that also draws the structure as an
  arc diagram, one color per stem.

## Citation

If you use fuzzyfold in your research, please cite:

> Stefan Badelt (2026): fuzzyfold: a high-performance framework for
> stochastic RNA folding kinetics. bioRxiv preprint,
> https://doi.org/10.64898/2026.06.17.732885
