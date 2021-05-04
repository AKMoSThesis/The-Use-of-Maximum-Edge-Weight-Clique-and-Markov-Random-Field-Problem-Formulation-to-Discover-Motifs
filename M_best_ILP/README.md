# M_best_ILP

M\_best\_ILP is a python package that finds the best motifs among multiple sequences using the ***ILP Formulation*** (based on [*A combinatorial optimization approach for diverse motif finding applications*](https://link.springer.com/article/10.1186/1748-7188-1-13)).

## Usage

1. Open terminal
2. Make sure you are in the ***/M\_best\_ILP/tests/*** directory
3. Into the console, type:

``` sh
$ python M_best_ILP_test.py [FILENAME] [MOTIF_LENGTH] [SOLVER] [SOLUTION_AMOUNT]
```

`[FILENAME]` = Sequence *filename* (sequences found in ***/M\_best\_ILP/seqs_and_data/***)

`[MOTIF_LENGTH]` = Desired *length* of motif

`[SOLVER]` = Desired solver (either *Gurobi* or *ORtools*)

`[SOLUTION_AMOUNT]` = If chosen `[SOLVER]` is Gurobi then specify the *number* of desired solutions, if chosen `[SOLVER]` is ORtools then specify whether multiple solutions are desired (*True*) or not (*False*)

A csv file with the name *M\_best\_ILP\_Solutions\_(SOLVER)\_(FILENAME).csv* will be generated in the ***/M\_best\_ILP/seqs\_and\_data/*** directory.
