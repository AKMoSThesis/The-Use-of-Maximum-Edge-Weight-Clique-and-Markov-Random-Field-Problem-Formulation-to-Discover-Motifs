#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:31:11 2021

@author: iamserious
"""

import sys
import os

#Grab Motif_finder and helpers
os.chdir("..")
os.chdir("../B_and_B_method/B_and_B_method")
sys.path.insert(0, os.path.abspath(os.curdir))

from important_functions import return_sequences
from important_functions import return_motifs
from important_functions import return_simple_score
from important_functions import return_background
from important_functions import return_score
############## B_and_B_method-specific ###############
from important_functions import return_graph_info
from important_functions import MECQ
from important_functions import Expand
from important_functions import Calc_seq_and_ub
from important_functions import N
from important_functions import N_thorough
#####################################################
from important_functions import return_motif_info
from important_functions import print_motif_solutions

from B_and_B_method import B_and_B_method

#Go to sequences
os.chdir("..")
os.chdir("../B_and_B_method/seqs_and_data")
