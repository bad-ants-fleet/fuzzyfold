# fuzzyfold - Python Interface

This repository provides (preliminary) Python bindings for **fuzzyfold's** RNA
energy evaluation and kinetic simulations.

The Python module exposes:

- Energy evaluation (only default RNA parameters)
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
For now, see the script: 
[example.py]

