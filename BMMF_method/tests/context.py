#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 15:23:58 2020

@author: iamserious
"""

import sys
import os

#Grab Motif_finder and helpers
os.chdir("..")
os.chdir("../BMMF_method/BMMF_method")
sys.path.insert(0, os.path.abspath(os.curdir))

from important_functions import return_sequences
from important_functions import return_motifs
from important_functions import return_simple_score
from important_functions import return_background
from important_functions import return_score
############ BMMF_method-specific #################
from important_functions import BMMmotif
from important_functions import BMMFmotif
##################################################
from important_functions import return_solution_score
from important_functions import return_motif_info
from important_functions import print_motif_solutions

from BMMF_method import BMMF_method

#Go to sequences
os.chdir("..")
os.chdir("../BMMF_method/seqs_and_data")
