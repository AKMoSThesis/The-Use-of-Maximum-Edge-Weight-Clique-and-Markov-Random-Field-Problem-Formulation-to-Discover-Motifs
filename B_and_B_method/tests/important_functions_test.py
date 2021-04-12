#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:30:22 2021

@author: iamserious
"""
from context import return_sequences
from context import return_motifs
from context import return_simple_score
from context import return_background
from context import return_score
############## INSERT B_and_B HERE ###############
from context import return_graph_info
from context import MECQ
from context import Expand
from context import Calc_seq_and_ub
from context import N
from context import N_thorough
##################################################
from context import return_motif_info
from context import print_motif_solutions

from B_and_B_method import B_and_B_method

#Test parameters
num_solutions = 2
multiple_solutions = True
node_count = 6
partition_count = 2
edge_count = 9
partition_nodes = [[0,3], [3,6]]
edges = [[0,3,1], [0,4,2], [0,5,3], [1,3,4], [1,4,5], [1,5,6], [2,3,7], [2,4,8], [2,5,9]]
weights = [1, 2, 3, 4, 5, 6, 7, 8, 9]

print(N(5, partition_nodes))
print()

C = {'nodes' : [], 'tot_weight' : 0}
S = [0,1,2,3,4,5]

PI, upper = Calc_seq_and_ub(C, S, node_count, partition_nodes, edges)
print(PI, upper)
print()

C = {'nodes' : [2], 'tot_weight' : 0}
S = [0,1,3,4,5]

PI, upper = Calc_seq_and_ub(C, S, node_count, partition_nodes, edges)
print(PI, upper)
print()

cliques = MECQ(node_count, partition_count, edge_count, partition_nodes, edges, weights, Cinit = {'nodes': [], 'tot_weight' : 0})
print(cliques)
print()

#Test parameters
num_solutions = 2
multiple_solutions = True
node_count = 6
partition_count = 3
edge_count = 12
partition_nodes = [[0,2], [2,4], [4,6]]
edges = [[0,2,9], [0,3,1], [0,4,9], [0,5,1], [1,2,1], [1,3,1], [1,4,1], [1,5,1], [2,4,9], [2,5,1], [3,4,1], [3,5,1]]
weights = [9,1,9,1,1,1,1,1,9,1,1,1]

cliques = MECQ(node_count, partition_count, edge_count, partition_nodes, edges, weights, Cinit = {'nodes': [], 'tot_weight' : 0})
print(cliques)
print()

#Test parameters
num_solutions = 2
multiple_solutions = True
node_count = 4
partition_count = 3
edge_count = 4
partition_nodes = [[0,1], [1,2], [2,4]]
edges = [[0,1,10],[0,2,10],[0,3,50],[1,2,10]]
weights = [10,10,50,10]

#NOTE this will not work correctly because you need to use N_thorough instead of N
cliques = MECQ(node_count, partition_count, edge_count, partition_nodes, edges, weights, Cinit = {'nodes': [], 'tot_weight' : 0})
print(cliques)
print()
