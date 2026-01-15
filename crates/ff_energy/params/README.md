# fuzzyfold's parameter files

These parameter files are taken from the **ViennaRNA** package and are used to
test consistency with ViennaRNA predictions.

Some files contain minor modifications compared to the original ViennaRNA
source files:

- **`rna_turner2004`**
  - Contains an additional `Misc` parameter.

- **General formatting changes**
  - `# END` statements were removed or replaced by `#END`.

- **`rna_turner2004_ext`**
  - A combined parameter set consisting of:
    - `rna_turner2004`
    - Additional special hairpin energies

