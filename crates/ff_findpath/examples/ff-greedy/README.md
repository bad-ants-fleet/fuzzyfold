# Greedy folding path heuristics and energy barrier calculation

## Overview

The program `ff-greedy` calculates a **folding trajectory** of
nucleic acid sequences from a starting secondary structure to a target secondary structure.
The program shows which base pair moves (deletion/insertion) are needed to fold from the start to the target structure given a greedy folding path heuristcs (Morgan-Higgs algorithm).
Each trajectory is shown as a sequence of structures with applied moves and associated
energies. The program then gives statistics about the energy barrier along this trajectory. 

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
ff-greedy -f short.txt
```

The program reads the input file `short.txt` with the command line argument `-f` and calculates the greedy folding path and corresponding energies. The output is shown on the command line.


---

## Example output

An example folding path output is shown below:

```
Folding Path found with greedy algorithm:
-----------------
AGCCAUGAGUGUAUAGUGGGCCUAU 	 applied move 	 energy
.(((..............))).... 	 del(0, 0)   -2.2 kcal/mol
.(((.(.........)..))).... 	 ins(5, 15)      -0.3 kcal/mol
.(((((.........)).))).... 	 ins(4, 16) 	 -1.9 kcal/mol
..((((.........)).))..... 	 del(1, 20) 	 1.4 kcal/mol
...(((.........)).)...... 	 del(2, 19) 	 4.3 kcal/mol
....((.........))........ 	 del(3, 18) 	 3.8 kcal/mol
...(((.........)))....... 	 ins(3, 17) 	 0.5 kcal/mol
..((((.........))))...... 	 ins(2, 18) 	 -2.7 kcal/mol
-----------------
Statistics:
Saddle energy: 4.30 kcal/mol, Barrier energy: 6.50 kcal/mol, Start energy: -2.20 kcal/mol, End energy: -2.70 kcal/mol
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

For better energy barrier estimations of folding paths see the example of the program
[`ff-findpath`](../ff-findpath/README.md) which makes use of the findpath algorithm.
