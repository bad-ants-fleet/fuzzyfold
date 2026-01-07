# Enumerate local minima, or other types of macrostates

## Input

Provide a sequence and a starting structure in the typical extednded FASTA format:

```fasta
>dld3
UCAGUCUUCGCUGCGCUGUAUCGAUUCGGUUUCAGUUUUUAUUGC
.((((....)))).((((........))))...............
```

*(The file name itself is arbitrary but indicates that this structure
corresponds to 'local minimum 3'.)*

```bash
cat dld3_lm3.fa | ff-locmin --delta 3
```

---

To familiarize yourself with other parameters, such as `--maxdist`, `--sorted`, etc.:

```bash
ff-timecourse --help
```


