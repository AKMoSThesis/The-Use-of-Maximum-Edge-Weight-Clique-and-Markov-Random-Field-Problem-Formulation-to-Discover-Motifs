#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:13:08 2020

@author: iamserious
"""

from context import return_sequences
from context import return_motifs
from context import return_simple_score
from context import return_background
from context import return_score
############# M_best_ILP-specific ###################
from context import return_graph_info
from context import Gurobi_ILPmotif
from context import ORtools_ILPmotif
#####################################################
from context import return_motif_info
from context import print_motif_solutions
from context import return_draw_info
from context import draw_graph

#Test parameters
num_solutions = 2
multiple_solutions = True
node_count = 6
partition_count = 2
edge_count = 9
partition_nodes = [[0,3], [3,6]]
edges = [[0,3,1], [0,4,2], [0,5,3], [1,3,4], [1,4,5], [1,5,6], [2,3,7], [2,4,8], [2,5,9]]
weights = [1, 2, 3, 4, 5, 6, 7, 8, 9]

#Test
ILPmotif_solutions = Gurobi_ILPmotif(num_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights)
print(ILPmotif_solutions)
ILPmotif_solutions = ORtools_ILPmotif(multiple_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights)
print(ILPmotif_solutions)
