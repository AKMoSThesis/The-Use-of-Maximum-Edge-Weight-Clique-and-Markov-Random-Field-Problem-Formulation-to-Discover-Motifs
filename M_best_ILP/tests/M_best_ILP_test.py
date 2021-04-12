#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:12:27 2020

@author: iamserious
"""

import sys

from context import Gurobi_M_best_ILP 
from context import ORtools_M_best_ILP 

#Test parameters
#################################
filename = sys.argv[1]
motif_length = int(sys.argv[2])
#################################
solver = sys.argv[3]
if solver == "Gurobi":
    num_solutions = int(sys.argv[4])
elif solver == "ORtools":
    multiple_solutions = sys.argv[4]
#################################
draw_graphs = False
#################################

##Example
#filename = "seq_hello_world.txt"
#motif_length = 15
#solver = "Gurobi"
#num_solutions = 13
#multiple_solutions = "True"
#draw_graphs = False

#Test Gurobi
if solver == "Gurobi":
    Gurobi_ILPmotif_solutions, motif_solutions = Gurobi_M_best_ILP(filename, motif_length, num_solutions, draw_graphs)
#    print(Gurobi_ILPmotif_solutions)
#    print(motif_solutions)
    
#Test ORtools
elif solver == "ORtools":
    ORtools_ILPmotif_solutions, motif_solutions = ORtools_M_best_ILP(filename, motif_length, multiple_solutions, draw_graphs)
#    print(ORtools_ILPmotif_solutions)
#    print(motif_solutions)

#print()
#print("Example of using graph dictionary: 0th variable found in the 4th solution = " + ILPmotif_solutions[4]["solution"][0])
#print("Example of using motif dictionary: motif selected from squirrel sequence in the 4th solution = " + motif_solutions[4]["solution"]["squirrel"])
#print()
