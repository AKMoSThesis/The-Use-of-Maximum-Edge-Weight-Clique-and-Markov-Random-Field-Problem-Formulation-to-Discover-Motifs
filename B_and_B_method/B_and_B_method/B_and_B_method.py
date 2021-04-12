#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:27:38 2021

@author: iamserious
"""

import time

from important_functions import return_sequences
from important_functions import return_motifs
from important_functions import return_background
############## B_and_B_method-specific #####################
from important_functions import return_graph_info
from important_functions import MECQ
############################################################
from important_functions import return_motif_info
from important_functions import print_motif_solutions

#Find motifs REMOVE BACKGROUND FROM PRINT_MOTIF_SOLUTIONS
def B_and_B_method(filename, motif_length, num_solutions):
    
    #Solution details for printing
    solution_details = [filename, motif_length, num_solutions]
    
    #START time calculation
    start_time = time.perf_counter()
    
    seq_name_list, seq_list = return_sequences(filename) 
    #print(seq_name_list)
    #print(seq_list)
    
    motifs = return_motifs(seq_list, motif_length)
    background = return_background(seq_list)
    #print(motifs)
    #print(background)
    
    node_count, partition_count, edge_count, partition_nodes, edges, weights, all_motifs = return_graph_info(motifs, background)
    #print(node_count)
    #print(partition_count)
    #print(edge_count)
    #print(partition_nodes)
    #print(edges)
    #print(weights)
    #print(all_motifs)
    
    cliques = MECQ(node_count, partition_count, edge_count, partition_nodes, edges, weights, num_solutions)
    #print(cliques)
    
    motif_solutions = return_motif_info(cliques, all_motifs, seq_name_list)
    #print(motif_solutions)
    
    #END time calculation
    end_time = time.perf_counter()
    total_time = end_time - start_time
    
    #Add to solution details
    solution_details.append(total_time)
    
    #Print motif solutions
    print_motif_solutions(motif_solutions, solution_details, background, print_location = "csv_file")
    
    return motif_solutions
