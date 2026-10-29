# Findpath folding path heuristics and energy barrier calculation

## Overview

The program `ff-findpath` calculates a **folding trajectory** of
nucleic acid sequences from a starting secondary structure to a target secondary structure.
The program shows which base pair moves (deletion/insertion) are needed to fold from the start to the target structure given a findpath folding path heuristcs which **improves** upon the greedy folding paths heuristics. 
Each trajectory is shown as a sequence of structures with applied moves and associated
energies. The program then gives statistics about the energy barrier along this trajectory. 

## Improving upon the Morgan-Higgs (greedy) algorithm

The algorithm can be seen as an extension of the idea coined by the Morgan-Higgs
heuristic. It introduces an additional parameter, the **search width** `m`, which
regulates the search by retaining all `m` best paths for every distance step. The
underlying principle is that optimal direct paths will almost never be the result
of purely greedy decisions, however, often the second or third best move from a
given intermediate may be sufficient for finding the best possible path. 

Practically, when the search width is set to `m = 1`, the resulting path is computed
greedily and yields the same outcome as `ff-greedy`.

---

## Using an input file

You can start a trajectory simulation from a predefined **input file**, which
contains a RNA sequence, a starting secondary structure and a target secondary structure (both given in a Dot-Bracket Notation).

For example, the file `short.txt` specifies a molecule in a particular
starting and target conformation:

```
AGCCAUGAGUGUAUAGUGGGCCUAU
.(((..............)))....
..((((.........))))......
```

You can run the greedy folding algorithm as follows:

```bash
ff-findpath -f short.txt -m 100
```

The program reads the input file `short.txt` with the command line arguments `-f` and a search width parameter `-m`and  calculates the findpath folding path and corresponding energies. The output is shown on the command line.


---

## Example output

An example folding path output is shown below:

```
Folding Path found with findpath algorithm:
-----------------
AGCCAUGAGUGUAUAGUGGGCCUAU 	 applied move 	 energy
.(((..............))).... 	 del(0, 0) 	 -2.2 kcal/mol
.((................)).... 	 del(3, 18) 	 1.6 kcal/mol
.(((.............).)).... 	 ins(3, 17) 	 1.4 kcal/mol
.((((...........)).)).... 	 ins(4, 16) 	 0.5 kcal/mol
.(((((.........))).)).... 	 ins(5, 15) 	 -1.9 kcal/mol
..((((.........))).)..... 	 del(1, 20) 	 1.4 kcal/mol
...(((.........)))....... 	 del(2, 19) 	 0.5 kcal/mol
..((((.........))))...... 	 ins(2, 18) 	 -2.7 kcal/mol
-----------------
Statistics:
Saddle energy: 1.60 kcal/mol, Barrier energy: 3.80 kcal/mol, Start energy: -2.20 kcal/mol, End energy: -2.70 kcal/mol```
```


Each line represents one **structure** visited during the trajectory, with columns:

| Header | Description |
|---------|-------------|
| **sequence** | Corresponding structure in dot-bracket notation. |
| **applied move** | Applied move to get to the corresponding structure. |
| **energy** | Free energy evaluation (kcal/mol). |

The summary statistics of the folding path are also shown:
The Saddle energy is the highest energy structure in the folding path. The Barrier energy is defined as the difference of the starting energy and the saddle energy.



---

## See also

For greedy energy barrier estimations of folding paths see the example of the program
[`ff-greedy`](../ff-greedy/README.md) which makes use of the Morgan-Higgs algorithm.
