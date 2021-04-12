#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:13:28 2020

@author: iamserious
"""

import sys
import os

#Grab Motif_finder and helpers
os.chdir("..")
os.chdir("../M_best_ILP/M_best_ILP")
sys.path.insert(0, os.path.abspath(os.curdir))

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

from M_best_ILP import Gurobi_M_best_ILP 
from M_best_ILP import ORtools_M_best_ILP 

#Go to sequences
os.chdir("..")
os.chdir("../M_best_ILP/seqs_and_data")
