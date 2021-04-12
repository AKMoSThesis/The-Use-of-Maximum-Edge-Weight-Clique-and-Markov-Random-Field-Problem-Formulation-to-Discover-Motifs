#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:02:25 2020

@author: iamserious
"""

import time

from important_functions import return_sequences
from important_functions import return_motifs
from important_functions import return_simple_score
from important_functions import return_background
from important_functions import return_score
############# M_best_ILP-specific ###################
from important_functions import return_graph_info
from important_functions import Gurobi_ILPmotif
from important_functions import ORtools_ILPmotif
#####################################################
from important_functions import return_motif_info
from important_functions import print_motif_solutions
from important_functions import return_draw_info
from important_functions import draw_graph

#Run Gurobi_M_best_ILP
def Gurobi_M_best_ILP(filename, motif_length, num_solutions, draw_graphs = False):
    
    #Solution details for printing
    solution_details = [filename, "Gurobi", motif_length, num_solutions]
    
    #START time calculation
    start_time = time.perf_counter()
    
    #Obtain all relevant info
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
    
    ILPmotif_solutions = Gurobi_ILPmotif(num_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights)
    #print(ILPmotif_solutions)
    
    motif_solutions = return_motif_info(ILPmotif_solutions, all_motifs, seq_name_list)
    #print(motif_solutions)
    
    #END time calculation
    end_time = time.perf_counter()
    total_time = end_time - start_time
    
    #Add to solution details
    solution_details.append(total_time)
    
    #Print motif solutions
    print_motif_solutions(motif_solutions, solution_details, print_location = "csv_file")
    
    #Draw graphs (if necessary)
    if draw_graphs == True:
        draw_solutions_nodes, draw_solutions_edges = return_draw_info(ILPmotif_solutions, edges)
        #print(draw_solutions_nodes)
        #print(draw_solutions_edges)
        draw_graph(draw_solutions_nodes, draw_solutions_edges, partition_count, edge_count, partition_nodes, edges)
    
    return ILPmotif_solutions, motif_solutions

#Run ORtools_M_best_ILP
def ORtools_M_best_ILP(filename, motif_length, multiple_solutions, draw_graphs = False):
    
    #Solution details for printing
    solution_details = [filename, "ORtools", motif_length]
    
    #START time calculation
    start_time = time.perf_counter()
    
    #Obtain all relevant info
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
    
    ILPmotif_solutions = ORtools_ILPmotif(multiple_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights)
    #print(ILPmotif_solutions)
    
    motif_solutions = return_motif_info(ILPmotif_solutions, all_motifs, seq_name_list)
    #print(motif_solutions)
    
    #Add to solution details
    solution_details.append(len(motif_solutions))
    
    #END time calculation
    end_time = time.perf_counter()
    total_time = end_time - start_time
    
    #Add to solution details
    solution_details.append(total_time)
    
    #Print motif solutions
    print_motif_solutions(motif_solutions, solution_details, print_location = "csv_file")
    
    #Draw graphs (if necessary)
    if draw_graphs == True:
        draw_solutions_nodes, draw_solutions_edges = return_draw_info(ILPmotif_solutions, edges)
        #print(draw_solutions_nodes)
        #print(draw_solutions_edges)
        draw_graph(draw_solutions_nodes, draw_solutions_edges, partition_count, edge_count, partition_nodes, edges)
    
    return ILPmotif_solutions, motif_solutions
