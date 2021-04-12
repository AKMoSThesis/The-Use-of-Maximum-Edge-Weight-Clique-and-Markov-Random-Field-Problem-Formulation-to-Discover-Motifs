7#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:29:05 2021

@author: iamserious
"""

import sys

from context import B_and_B_method

#Test parameters
#################################
filename = sys.argv[1]
motif_length = int(sys.argv[2])
#################################
num_solutions = int(sys.argv[3])
#################################

##Example
#filename = "seq_hello_world.txt"
#motif_length = 15
#num_solutions = 13

#Test B_and_B
motif_solutions = B_and_B_method(filename, motif_length, num_solutions)
#print(motif_solutions)
