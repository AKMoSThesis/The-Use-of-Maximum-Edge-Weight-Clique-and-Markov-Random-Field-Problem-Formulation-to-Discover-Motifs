#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 13:20:20 2020

@author: iamserious
"""

import numpy as np
import math
from loopy_BP import *
import csv

#Return list of sequences
def return_sequences(filename):

    #Read raw file
    file = open(filename, 'r')
    file_info = file.readlines()
    file.close()

    #Get seq dict and seq list from file
    seq_name_list = []
    seq_list = []
    while file_info[0][0] == '>':
        seq_name = file_info[0].replace(">", "")
        seq_name = seq_name.replace("\n", "")
        file_info.pop(0) #Remove the seq name to move on
        seq = ""
        while file_info[0][0] != '>' and len(file_info) > 1:
            seq += file_info[0].replace("\n", "") #Add seq group to seq
            file_info.pop(0) #Remove seq group to move on
        if len(file_info) == 1:
            seq += file_info[0] #When last line reached, add final seq group to seq
        seq = seq.upper() #Capitalize all letters
        seq_name_list.append(seq_name) #Add seq name to seq name list
        seq_list.append(seq) #Add seq to seq list

    return seq_name_list, seq_list

#Return list of motifs
def return_motifs(seq_list, motif_length):

    motifs = []
    for seq in seq_list:
        seq_motifs = []
        for num in range(len(seq) - motif_length + 1):
            seq_motifs.append(seq[num : num+motif_length]) #Add seq motif to list of seq motifs
        motifs.append(seq_motifs) #Add list of seq motifs to list of all motifs

    return motifs

#Return simple score of two motifs
def return_simple_score(motif1, motif2):

    length = len(motif1) #Find length of motifs
    score = 0
    for n in range(length): #Look through all aligned bases of both motifs
        if motif1[n] == motif2[n]: #Add to score iff the bases align
            score += 1

    return score

#Return background of sequences
def return_background(seq_list):

    background = {"total" : 0}
    for seq in seq_list:
        for base in seq: #Look through every existing base of all sequences
            if base not in background: #Add base to background if it doesn't exist
                background[base] = 0
            background[base] += 1
            background["total"] +=1

    return background

#Return score of two motifs
def return_score(motif1, motif2, background):

    length = len(motif1) #Find length of motifs
    score = 0
    for n in range(length): #Look through all aligned bases of both motifs
        if motif1[n] == motif2[n]: #Add to score iff the bases align
            f = background[motif1[n]] / background["total"]
            score += math.log(1/(f))

    return score

##################################################### BMMF_method-specific ################################################################

#Find Best MM of motifs
def BMMmotif(motifs, background, constraints = []):

    num_sequences = len(motifs)

    #Considering the constraints
    motifs_constr = []
    for seq_motifs in motifs:
        seq_motifs_constr = []
        for motif_constr in seq_motifs:
            seq_motifs_constr.append(motif_constr)
        motifs_constr.append(seq_motifs_constr)
    for constr in constraints:
        i_constr, j_constr, equality = constr.values()
        if equality:
            motifs_constr[i_constr] = [j_constr]
        else:
            motifs_constr[i_constr].remove(j_constr)

    #Create mrf
    mrf = factor_graph()

    #Create factors for mrf
    factor_num = 0
    for i1 in range(num_sequences):
        for i2 in range(i1 + 1, num_sequences):

            #Making the factor arrays
            factors_list = []
            for motif1 in motifs_constr[i1]:
                factor_list = []
                for motif2 in motifs_constr[i2]:
#                    score = return_simple_score(motif1, motif2)
                    score = return_score(motif1, motif2, background)
                    factor_list.append(score)
                factors_list.append(factor_list)

            #Making the factors
            factor_motif = factor(['x' + str(i1), 'x' + str(i2)], np.array(factors_list))
            mrf.add_factor_node('f' + str(factor_num), factor_motif)
            factor_num += 1

    #Create loopy belief
    lbp = loopy_belief_propagation(mrf)

    #Finding the best motifs
    motif_likelyhood = 0
    BMMmotifs = []
    for i in range(num_sequences):
        motif_belief = lbp.belief('x' + str(i), 10).get_distribution()
        motif_index = np.argmax(motif_belief)
        motif_likelyhood += np.max(motif_belief)
        BMMmotifs.append(motifs_constr[i][motif_index])

    return BMMmotifs, motif_likelyhood

#Perform BMMF
def BMMFmotif(motifs, background, num_solutions):

    #Initialize for m (list of solutions)
    l = len(motifs) #Motif length
    m = [] #Solutions list
    constraints = [] #Constraints for each solution

    #Find all mts
    for t in range(num_solutions): #Number of solutions

        solution_dict = {"solution" : [], "score" : 0}

        #Add first solution
        if t == 0:

            BMMmotifs, motif_likelyhood = BMMmotif(motifs, background)
            solution_dict["solution"] = BMMmotifs
            solution_dict["score"] = motif_likelyhood
            m.append(solution_dict)
            constraints.append([])

        #Add all solutions
        else:

            #Find the max i,j,s
            constraint_likelyhood = 0
            max_vars = (-1, "", -1)
            for s in range(len(m)):
                for i in range(l):
                    s_constraints = []
                    for s_constraint in constraints[s]:
                        s_constraints.append(s_constraint)
                    s_constraints.append({"i" : i, "j" : m[s]["solution"][i], "equality" : False})
                    #Check if constraints repeat
                    skip = False
                    for c1 in range(len(s_constraints)):
                        for c2 in range(c1 + 1, len(s_constraints)):
                            if (s_constraints[c1]["i"] == s_constraints[c2]["i"]) and (s_constraints[c1]["j"] == s_constraints[c2]["j"]) and (s_constraints[c1]["equality"] != s_constraints[c2]["equality"]):
                               skip = True
                    if skip:
                        pass
                    else:
                        BMMmotifs, motif_likelyhood = BMMmotif(motifs, background, s_constraints)
                        if motif_likelyhood > constraint_likelyhood:
                            max_vars = (i, BMMmotifs[i], s)
                            constraint_likelyhood = motif_likelyhood

            #Edit constraints
            i, j, s = max_vars
            t_constraints = []
            for s_constraint in constraints[s]:
                t_constraints.append(s_constraint)
            t_constraints.append({"i" : i, "j" : j, "equality" : True})
            constraints.append(t_constraints)
            constraints[s].append({"i" : i, "j" : j, "equality" : False})

            #Find the next solution
            BMMmotifs, motif_likelyhood = BMMmotif(motifs, background, constraints[t])
            solution_dict["solution"] = BMMmotifs
            solution_dict["score"] = motif_likelyhood
            m.append(solution_dict)

    return m

###########################################################################################################################################

#Calculates solution score
def return_solution_score(motif_solution, background):

    motifs = []
    for motif in motif_solution["solution"].values():
        motifs.append(motif)

    solution_score = 0
    for i in range(len(motifs)):
        for j in range(i + 1, len(motifs)):
            solution_score += return_score(motifs[i], motifs[j], background)

    return solution_score

#Return motif info from BMMF solution
def return_motif_info(m, seq_name_list):

    motif_solutions = []
    for m_sol in m:
        motif_solution = {}
        seq_count = 0
        for motif in m_sol["solution"]: #Look through all solutions of BMMF solution
            motif_solution[seq_name_list[seq_count]] = motif #Add selected motif to dict of one motif solution
            seq_count += 1
        score = m_sol["score"] #Return the score of the solution
        motif_solutions.append({"solution" : motif_solution, "score" : score})

    return motif_solutions

#Print the motif solutions in console (REMOVE BACKGROUND INPUT, AND REMOVE RETURN_SOLUTION_SCORE)
def print_motif_solutions(motif_solutions, solution_details, background, print_location = "console", print_format = "elegant"):

    #Unpack solution details
    filename, motif_length, num_solutions, total_time = solution_details

    #Print the motif solutions in the console
    if print_location == "console":

        print("\nSOLUTIONS GENERATED!")
        print("filename =", filename)
        print("motif_length =", motif_length)
        print("num_solutions =", num_solutions)
        print("total_time =", total_time)

        #Print the motif solutions in an easy-to-read-format
        if print_format == "elegant":

            print("-" * 60)
            for n in range(len(motif_solutions)):
                print("-" * 60)
                print("Solution #" + str(n) + " choses the following motifs:")
                print("-" * 60 + "\n")
                print("{:<30s}{:>10s}".format("SEQUENCE", "SELECTED MOTIFS\n"))
                for seq_name, motif in motif_solutions[n]["solution"].items():
                    print("{:<30s}{:>10s}".format(seq_name, motif))
#                print("\n\nScore of solution (LOOPY BP SCORE) = " + str(motif_solutions[n]["score"]) + "\n")
                print("\n\nScore of solution = " + str(return_solution_score(motif_solutions[n], background)) + "\n") 
            print("-" * 60)
            print("-" * 60)
            print("")

        #Print the motif solutions in an easy-to-record-format
        elif print_format == "testing":

            print("-" * 60)
            for n in range(len(motif_solutions)):
                print("-" * 60)
                print("Solution #" + str(n) + " choses the following motifs:")
                print("-" * 60 + "\n")
                print("SELECTED MOTIFS\n")
                for seq_name, motif in motif_solutions[n]["solution"].items():
                    print(motif)
                print(return_solution_score(motif_solutions[n], background), "\n")
            print("-" * 60)
            print("-" * 60)
            print("")

    #Print the motif solutions in a text file
    elif print_location == "txt_file":

        solution_filename = "BMMF_Solutions_(" + filename.replace(".txt", "") + ")" + ".txt"
        motif_file = open(solution_filename, "a")
        print("SOLUTIONS GENERATED!", file = motif_file)
        print("motif_length =", motif_length, file = motif_file)
        print("num_solutions =", num_solutions, file = motif_file)
        print("total_time =", total_time, file = motif_file)

        #Print the motif solutions in an easy-to-read-format
        if print_format == "elegant":

            print("-" * 60, file = motif_file)
            for n in range(len(motif_solutions)):
                print("-" * 60, file =  motif_file)
                print("Solution #" + str(n) + " choses the following motifs:", file = motif_file)
                print("-" * 60 + "\n", file =  motif_file)
                print("{:<30s}{:>10s}".format("SEQUENCE", "SELECTED MOTIFS\n"), file = motif_file)
                for seq_name, motif in motif_solutions[n]["solution"].items():
                    print("{:<30s}{:>10s}".format(seq_name, motif), file = motif_file)
                print("\n\nScore of solution = " + str(return_solution_score(motif_solutions[n], background)) + "\n", file = motif_file)
            print("-" * 60, file = motif_file)
            print("-" * 60, file = motif_file)
            print("\n\n\n\n\n", file = motif_file)

        #Print the motif solutions in an easy-to-record-format
        elif print_format == "testing":

            print("-" * 60, file = motif_file)
            for n in range(len(motif_solutions)):
                print("-" * 60, file = motif_file)
                print("Solution #" + str(n) + " choses the following motifs:", file = motif_file)
                print("-" * 60 + "\n", file = motif_file)
                print("SELECTED MOTIFS\n", file =  motif_file)
                for seq_name, motif in motif_solutions[n]["solution"].items():
                    print(motif, file = motif_file)
                print(return_solution_score(motif_solutions[n], background), "\n", file = motif_file)
            print("-" * 60, file = motif_file)
            print("-" * 60, file = motif_file)
            print("\n\n\n\n\n", file = motif_file)

        motif_file.close()

    #Print the motif solutions in a csv file
    elif print_location == "csv_file":

        solution_filename = "BMMF_Solutions_(" + filename.replace(".txt", "") + ")" + ".csv"
        motif_file = open(solution_filename, "a")
        motif_writer = csv.writer(motif_file)
        motif_writer.writerow(["Motif length:", motif_length])
        motif_writer.writerow(["num_solutions:", num_solutions])
        motif_writer.writerow(["Time needed:", total_time])

        header = ["Motifs"]
        for n in range(len(motif_solutions)):
            header.append("Solution " + str(n))
        motif_writer.writerow(header)

        for seq_name in motif_solutions[0]["solution"].keys():
            seq_row = [seq_name]
            for motif_solution in motif_solutions:
                seq_row.append(motif_solution["solution"][seq_name])
            motif_writer.writerow(seq_row)

        score_row = ["SCORE"]
        for motif_solution in motif_solutions:
            score_row.append(return_solution_score(motif_solution, background))
        motif_writer.writerow(score_row)

        motif_writer.writerow([])
        motif_writer.writerow([])

        motif_file.close()
