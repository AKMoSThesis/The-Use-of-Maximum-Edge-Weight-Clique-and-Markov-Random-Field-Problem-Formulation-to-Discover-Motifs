# B_and_B_method

B\_and\_B\_method is a python package that finds the best motifs among multiple sequences using the ***branch-and-bound algorithm*** (based off of [*A maximum edge-weight clique extraction algorithm based on branch-and-bound*](https://www-sciencedirect-com.proxy.library.carleton.ca/science/article/pii/S1572528620300177?via%3Dihub)).

## Usage

1. Open terminal
2. Make sure you are in the ***/B\_and\_B\_method/tests/*** directory
3. Into the console, type:

``` sh
$ python B_and_B_method_test.py [FILENAME] [MOTIF_LENGTH] [SOLUTION_AMOUNT]
```

`[FILENAME]` = Sequence *filename* (sequences found in ***/B\_and\_B\_method/seqs_and_data/***)

`[MOTIF_LENGTH]` = Desired *length* of motif

`[SOLUTION_AMOUNT]` = Specify the *number* of desired solutions

A csv file with the name *B\_and\_B_method\_Solutions\_(FILENAME).csv* will be generated in the ***/B\_and\_B\_method/seqs\_and\_data/*** directory.
