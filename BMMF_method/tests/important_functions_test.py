#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:23:06 2020

@author: iamserious
"""

from context import return_sequences
from context import return_motifs
from context import return_simple_score
from context import return_background
from context import return_score
############ BMMF_method-specific #################
from context import BMMmotif
from context import BMMFmotif
##################################################
from context import return_solution_score
from context import return_motif_info
from context import print_motif_solutions

filename = "seq_hello_world.txt"
motif_length = 15
i_const = 2
j_const = "AATGAATGCGATCGA"
j_const_z = "ZZZZZZZZZZZZZZZ"
iterations = 3
a_constraint = {"i" : i_const, "j" : j_const, "equality" : True}
b_constraint = {"i" : i_const, "j" : j_const_z, "equality" : False}

seq_name_list, seq_list = return_sequences(filename) 
print(seq_name_list)
print(seq_list)
    
motifs = return_motifs(seq_list, motif_length)
print(motifs)

background = return_background(seq_list)
print(background)

BMMmotifs, motif_likelyhood = BMMmotif(motifs, background)
print(BMMmotifs)
print(motif_likelyhood)

BMMmotifs, motif_likelyhood = BMMmotif(motifs, background, [a_constraint])
print(BMMmotifs)
print(motif_likelyhood)

BMMmotifs, motif_likelyhood = BMMmotif(motifs, background, [b_constraint])
print(BMMmotifs)
print(motif_likelyhood)

#print("\nChecking BMMF!!!\n")
#m = BMMFmotif(motifs, background, 0)
#print("m =", m)

#print("\nChecking BMMF!!!\n")
#m = BMMFmotif(motifs, background, 1)
#print("m =", m)

print("\nChecking BMMF!!!\n")
m = BMMFmotif(motifs, background, 3)
print("m =", m)
