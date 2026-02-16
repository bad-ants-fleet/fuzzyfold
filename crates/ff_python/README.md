# fuzzyfold - Python Interface (draft)

This repository provides Python bindings for the Rust-based **fuzzyfold**
nucleic acid folding and kinetic simulation engine.

The computational core is implemented in Rust.
The Python layer is a thin wrapper around the Rust engine.

The Python module exposes:

- Energy evaluation (`ViennaRNA`)
- Stochastic kinetic simulation (`Simulator`)

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



