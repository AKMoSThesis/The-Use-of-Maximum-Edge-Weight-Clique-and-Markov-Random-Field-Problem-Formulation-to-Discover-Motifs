#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 10:36:20 2020

@author: iamserious
"""

import time

from important_functions import return_sequences
from important_functions import return_motifs
from important_functions import return_background
############# BMMF_method-specific ###################
from important_functions import BMMFmotif
######################################################
from important_functions import return_motif_info
from important_functions import print_motif_solutions

#Run BMMF_method 
def BMMF_method(filename, motif_length, num_solutions):
    
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
    
    m = BMMFmotif(motifs, background, num_solutions)
    #print(m)
    
    motif_solutions = return_motif_info(m, seq_name_list)
    #print(motif_solutions)
    
    #END time calculation
    end_time = time.perf_counter()
    total_time = end_time - start_time
    
    #Add to solution details
    solution_details.append(total_time)
    
    #Print motif solutions
    print_motif_solutions(motif_solutions, solution_details, background, print_location = "csv_file")
    
    return m, motif_solutions
