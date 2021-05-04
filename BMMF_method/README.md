# BMMF_method

BMMF\_Method is a python package that finds the best motifs among multiple sequences using the ***loopy belief propagation*** (based on [*M are better than one: an ensemble-based motif finder and its application to regulatory element prediction*](https://academic.oup.com/bioinformatics/article/25/7/868/211358)).

## Usage

1. Open terminal
2. Make sure you are in the ***/BMMF\_method/tests/*** directory
3. Into the console, type:

``` sh
$ python BMMF_method_test.py [FILENAME] [MOTIF_LENGTH] [SOLUTION_AMOUNT]
```

`[FILENAME]` = Sequence *filename* (sequences found in ***/BMMF\_method/seqs\_and\_data/***)

`[MOTIF_LENGTH]` = Desired *length* of motif

`[SOLUTION_AMOUNT]` = Specify the *number* of desired solutions

A csv file with the name *BMMF\_method\_Solutions\_(FILENAME).csv* will be generated in the ***/BMMF\_method/seqs\_and\_data/*** directory.
