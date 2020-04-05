# proof_graph_analyzer

Scripts to analyze proof graphs and perform mergeability computations.

## Community File Inputs
For all community (cmty) file inputs, the file should have one var (zero-based) and its community on each line, e.g.,
```
0 1
1 2
2 1
```
The above example indicates that variable 0 (which corresponds to 1 in the dimacs) is in community 1, variable 1 (aka 2 in dimacs) is in community 2, and so on.

Louvain scripts in other supplemental `backdoors` repos should facilitate this.

## merge_generator.cpp
```
USAGE: ./merge_generator cnf_file out_cnf_file intermediate_flips max_flips reduce0_increase1
intermediate_flips (N) -- after every N literal flips that occur, dump the current CNF file. This allows us to take snapshots of the CNF as its mergeability increases (or decreases).
max_flips (N) -- stop flipping literals after N flips and dump the final CNF
reduce0_increase1 (0 or 1) -- if set to 0, reduce mergeability, else increase.

```

 Takes a cnf and greedily increases merges in the formula.
 
 Defintions:
 *   enabler literal - for a clause pair, resolution over the literal allows a merge to happen over some other literal
 *   merge literal - there exists a clause pair that has a merge over the literal
 *   double resolution - a clause pair resolves on two literals, and merges on a third: (a v b v c) ^ (!a v !b v c)

 The program takes an input CNF. The clause "skeleton" is "locked", i.e., the set of clauses remain the same,
 but the polarities of literals can be changed. The number of each literal is also locked, so in order to
 swap v to !v, !v must be swapped to v elsewhere. This ensures that the community structure and literal popularity
 distribution of the original instance is maintained. It also ensures that if 2 clauses are already merge+res, they
 will still be.

 If a literal is an enabler or a merge, it is also locked. If there is a double resolution, one may flip.
 
## merge_resolution_dimacs_check.cpp

```
USAGE: ./proof_graph_analyzer cnf_file cmty_file out_file
```

Outputs mergeability data for a CNF file such as:
```
ClausePairs
NumResolutions
NumMerges
NumMergePairs
AvgCmtyMergeOverResolutions
NormalizedMerges
NormalizedResolutions
```
