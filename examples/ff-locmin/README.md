# Enumerate local minima, or other types neighboring structures

## Example

The file `dld3_lm3.fa` provides a sequence and a starting structure in the following format:

```fasta
>lm3 (dld3)
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
.((((....)))).((((........))))...............
```

Use this file as input to enumerate the neighborhood using base-pair moves. The following
command returns all reachable structures with energy smaller or equal E(s) + 2.80 kcal/mol:

```bash
cat dld3_lm3.fa | ff-locmin --delta 2.8 --sorted
```


```
>lm3 (dld3) (delta = 2.80)
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
.((((....)))).((((........))))............... -3.90
.(((......))).((((........))))............... -2.60
.((((....))))..(((........)))................ -2.40
.((((....)))).((((.(.....)))))............... -1.90
..(((....)))..((((........))))............... -1.60
.(((......)))..(((........)))................ -1.10
.(((.(...)))).((((........))))............... -1.10
.((((....)))).((((.(....).))))............... -1.10
```


---

Be careful, the energy delta is relative to the starting structure. If it is
too high, then many more low-energy structures will be enumerated. Use `--help`
to familiarize yourself with the other parameters, such as maximum base-pair
distances.

```bash
ff-timecourse --help
```


