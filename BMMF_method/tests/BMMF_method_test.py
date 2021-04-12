#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 13:13:05 2020

@author: iamserious
"""

import sys

from context import BMMF_method

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

#Test BMMF
m, motif_solutions = BMMF_method(filename, motif_length, num_solutions)
#print(m)
#print(motif_solutions)
