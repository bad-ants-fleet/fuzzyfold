# Single RNA folding trajectories from stochastic simulations

## Overview

The program `ff-trajectory` simulates **individual folding trajectories** of
nucleic acid sequences, showing how secondary structures evolve over time.
Each trajectory is represented as a sequence of structures with associated
energies and transition times.

---

## Using FASTA input

You can start a trajectory simulation from a predefined **FASTA file**, which
contains a sequence and optionally an initial structure.

For example, the file `dld3.in` defines an open-chain configuration:

```fasta
>dld3
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
```

while the file `dld3_lm3.in` specifies the same sequence but with a folded
starting conformation:

```fasta
>dld3
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
.((((....)))).((((........))))...............
```

You can run the trajectory simulation for either case as follows:

```bash
cat dld3.in | ff-trajectory --t-end 50
```

or

```bash
cat dld3_lm3.in | ff-trajectory --t-end 50
```

The simulation continues until **time = 50** arbitrary time units. (See rate
model details for conversion to wall-clock time.) All visited structures are
printed with their energies and transition times.

---

## Quick example with a random sequence

For a quick test, you can generate a random RNA sequence and simulate its trajectory in one command:

```bash
ff-randseq -l 50 | ff-trajectory --t-end 50
```

or, when running from the repository:

```bash
cargo run --bin ff-randseq -- -l 50 | cargo run --bin ff-trajectory -- --t-end 50
```

This creates a random 50-nt sequence and follows its stochastic folding trajectory until 50 seconds.

---

## Example output

An example asimulation output is shown below:

```
GCGUUUCCAGGGUUUAGACGGACGGGUGUGACUCGCCCAGCCCCGACCUC   energy   arrival time   waiting time    mean-waiting
..................................................     0.00   0.00000000e0  5.48838308e-1    1.15990269e0
.......................(.......)..................     2.10  5.48838308e-1  7.42186429e-1   6.00461698e-1
.......................(.......).....(......).....     4.50   1.29102474e0  6.62519714e-1   3.57034800e-1
.....................................(......).....     2.40   1.95354445e0  8.70308936e-3   5.15932279e-1
..................................................     0.00   1.96224754e0   2.24827455e0    1.15990269e0
..................(......)........................     1.90   4.21052209e0   1.08942365e0   4.14255646e-1
.................((......)).......................     0.10   5.29994574e0  3.93454939e-1   5.63963354e-1
...............(.((......)).).....................    -0.20   5.69340068e0  1.06769149e-1   4.52730587e-1
.............(.(.((......)).).)...................     1.80   5.80016983e0  4.33243095e-1   3.16296938e-1
.............(((.((......)).)))...................    -2.40   6.23341293e0   1.28813977e0   8.21431625e-1
...........(.(((.((......)).))))..................    -1.60   7.52155270e0  2.45415452e-1   2.40731917e-1
.........(.(.(((.((......)).)))).)................    -1.50   7.76696815e0  6.94202383e-2   4.10851300e-1
.........(((.(((.((......)).))))))................    -5.60   7.83638839e0   3.47571920e0    3.26611240e0
.........(((..((.((......)).)).)))................    -2.90   1.13121076e1   1.55874453e0   4.13337266e-1
.........(((.(((.((......)).))))))................    -5.60   1.28708521e1   2.36328447e0    3.26611240e0
.........(((.(((.((......)).))))))(.....).........    -2.10   1.52341366e1  7.41356896e-3   4.84982503e-1
.........(((.(((.((......)).))))))((...)).........    -4.70   1.52415502e1   2.00829892e1    1.23961105e1
..(......(((.(((.((......)).))))))((...))......)..    -0.50   3.53245394e1   1.26546682e0   7.99897165e-1
.........(((.(((.((......)).))))))((...)).........    -4.70   3.65900062e1   4.58067121e1    1.23961105e1
```

Each line represents one **structure** visited during the trajectory, with columns:

| Column | Description |
|---------|-------------|
| **structure** | Secondary structure in dot-bracket notation |
| **energy** | Free energy (kcal/mol) |
| **arrival time** | Absolute time when the structure was reached |
| **waiting time** | Duration spent in that structure before transition |
| **mean-waiting** | Mean waiting time in this structure (1/flux) |

---

## See also

For ensemble-level analysis across many trajectories, see the example in
[`ff-timecourse`](../ff-timecourse/README.md), which aggregates population data
over multiple stochastic runs.


