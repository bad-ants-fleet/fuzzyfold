# fuzzyfold - Python Interface (draft)

This repository provides Python bindings for the Rust-based **fuzzyfold**
nucleic acid folding and kinetic simulation engine.

The computational core is implemented in Rust.
The Python layer is a thin wrapper around the Rust engine.

The Python module exposes:

- Energy evaluation (`ViennaRNA`)
- Stochastic kinetic simulation (`Simulator`)
- Pseudoknot parsing and energies (`region_tree`, `pseudo_energy`)

## Installation (Development)

The Python extension is built using `maturin`.

### 1. Install maturin

```bash
pip install maturin
```

### 2. Build and install locally
From the Python binding crate directory:

```bash
maturin develop --release
```


### 3. Basic usage

```python
import fuzzyfold as ff

seq = "UGCCUAGAGAGUCAGGUGAU"
db1 = ".((((.((...))))))..."

# Energy model
emodel = ff.ViennaRNA()
print(emodel.energy_of_structure(seq, db1))

# Kinetic simulator
ssa = ff.Simulator(k0=1)

for line in ssa.simulate(
        seq,
        None,
        t_ext=4000,
        t_end=4000,
        silent=False):
    print(line)
```


## Recommended alternative: uv

[uv](https://docs.astral.sh/uv/) replaces steps 1 and 2 above. It resolves and
locks quickly, and reads the `pyproject.toml` beside `Cargo.toml`:

```bash
uv venv --python 3.12          # create the environment
source .venv/bin/activate
uv pip install maturin         # the build backend
maturin develop --release      # run from this crate directory
```

`maturin develop` installs into the *active* environment, hence the `activate`
step. Having `pyproject.toml` present is also what makes the installed
distribution `fuzzyfold` matching the import name.

## Pseudoknots

Pseudoknotted structures use extended dot-bracket notation (`([{<` and `ABCD`).
Parsing follows Rastegari & Condon; energies use the Dirks-Pierce model.

| Export | Purpose |
|---|---|
| `region_tree(structure)` | `(pair_table, region_nodes)` — the closed-region tree |
| `region_tree_steps(structure)` | the same parse, emitted step by step |
| `pseudo_energy(seq, structure)` | free energy in kcal/mol |
| `pseudo_energy_loops(seq, structure)` | per-loop energy breakdown |
| `RegionNode`, `TreeStep` | the types those functions return |

`region_tree_steps` replays the construction rather than only reporting its
result, which is what drives the notebook visualisations.

See `example.py` for simple simulations and the energy breakdown of individual loops of an RNA pseudoknot.
