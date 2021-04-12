#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 20 09:28:07 2021

@author: iamserious
"""

import numpy as np
import math
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

##################################################### B_and_B_method-specific ################################################################

#Return graph info from motifs
def return_graph_info(motifs, background):
    
    #Collect node_count,  partition_count, partition_nodes, all_motifs
    node_count = 0 
    partition_count = 0
    partition_nodes = []
    all_motifs = []
    for seq_motifs in motifs: 
        first_node = node_count #First node of partition
        for motif in seq_motifs: #Look at each existing motif
            node_count += 1
            all_motifs.append(motif) #Add motif to list of all motifs
        last_node = node_count #Last node of partition
        partition_count += 1
        partition_nodes.append([first_node, last_node]) #Specify boundary nodes of partition
    
    #Collect edge_count, edges=list([first node, second_node, weight]), weights
    edge_count = 0
    edges = []
    weights = []
    for i in range(partition_count): 
        for j in range(i + 1, partition_count): #Look at all i,j partition combos
            first_node_i, last_node_i = partition_nodes[i] 
            first_node_j, last_node_j = partition_nodes[j]
            for u in range(first_node_i, last_node_i): #Look at all nodes in partition i
                for v in range(first_node_j, last_node_j): #Look at all nodes in partition j
                    motif1, motif2 = all_motifs[u], all_motifs[v]
                    weight = return_score(motif1, motif2, background) #Calculate weight of edge
                    edge_count += 1
                    edges.append([u, v, weight]) #Add edge=[first node, second_node, weight] to list of edges
                    weights.append(weight) #Add weight to list of weights
    
    return node_count, partition_count, edge_count, partition_nodes, edges, weights, all_motifs

#MECQ
def MECQ(node_count, partition_count, edge_count, partition_nodes, edges, weights, num_solutions, Cinit = {'nodes': [], 'tot_weight' : 0}):
    
    #Create Cmax
    Cmax = Cinit
    #Create S that contains all nodes that share an edge with all nodes in Cmax, not including the nodes themselves (note: we DON'T need to specify this in the loop)
    S = []
    for snode in range(node_count): #Look through every existing node
        neighbor_count = 0
        for cnode in Cmax['nodes']: 
            if snode in N(cnode, partition_nodes): #Check if the node is a neighbor to one of the nodes in Cmax
                neighbor_count += 1
        if neighbor_count == len(Cmax['nodes']): #Make sure that the node is a neighbor to all the nodes in Cmax (if it is, then neighbor_count should equal the number of nodes in Cmax)
            S.append(snode)    
    
    #Update Cmax, 
    Cmax = Expand(Cinit, Cmax, S, node_count, partition_nodes, edges)  
    
    #Find the number of desired solutions 
    cliques = [Cmax]
    for num in range(num_solutions-1):
        Cmax = Cinit
        Cmax = Expand(Cinit, Cmax, S, node_count, partition_nodes, edges, cliques)
        cliques.append(Cmax)
            
    return cliques 

#Expand
def Expand(C, Cmax, S, node_count, partition_nodes, edges, cliques = []): 
    
    #Return the superior Cmax when there are no longer any points left to add to the clique
    if S == []:
        if C['tot_weight'] > Cmax['tot_weight'] and C not in cliques: #This checks whether the C found is better than the current Cmax, and whether it is a clique that has been found before
            Cmax = {'nodes' : [], 'tot_weight' : 0}
            for node in C['nodes']:
                Cmax['nodes'].append(node)
            Cmax['tot_weight'] = C['tot_weight']
        return Cmax
    
    #Return PI (the node sequence) and upper (the upper bound of each node in the sequence) using node coloring to find a list of nodes to add
    PI, upper = Calc_seq_and_ub(C, S, node_count, partition_nodes, edges)
    
    #Find the next best node to add to the clique
    for p in PI: #Check through all desireable nodes to add to the clique
        if C['tot_weight'] + upper[p] >= Cmax['tot_weight']: #Check if adding the node increases the weight of the clique
            
            #Create new clique that has the new node added to the old clique
            Cnew = {'nodes' : [], 'tot_weight' : 0}
            Cnew['tot_weight'] = C['tot_weight'] 
            for node in C['nodes']:
                Cnew['nodes'].append(node) #Add every node in the old clique to the new clique
                for edge in edges: 
                    if node < p:
                        if [node, p] == edge[0:2]:
                            Cnew['tot_weight'] += edge[2] #Add the weight of the edge that connects the new node to each old node
                    elif p < node:
                        if [p, node] == edge[0:2]:
                            Cnew['tot_weight'] += edge[2] #Add the weight of the edge that connects the new node to each old node
            Cnew['nodes'].append(p) #Add the new node to the new clique
            
            #Create a new S that contains all nodes that share an edge with all nodes in Cnew, not including the nodes themselves (note: we DON'T need to specify this in the loop)
            Snew = []
            for snode in S: #Look through every node already in S
                neighbor_count = 0
                for cnode in Cnew['nodes']:
                    if snode in N(cnode, partition_nodes): #Check if the node is a neighbor to one of the nodes in Cmax
                        neighbor_count += 1
                if neighbor_count == len(Cnew['nodes']): #Make sure that the node is a neighbor to all the nodes in Cmax (if it is, then neighbor_count should equal the number of nodes in Cmax)
                    Snew.append(snode)            
        
            #Recurse Expand
            Cmax = Expand(Cnew, Cmax, Snew, node_count, partition_nodes, edges, cliques)
            
    return Cmax

#Calc_seq_and_ub
def Calc_seq_and_ub(C, S, node_count, partition_nodes, edges):
    
    #Create sigma (a list of weights)
    sigma = []
    for node in range(node_count):
        if node not in S:
            sigma.append(None) #If the node is already part of the clique, it will never also be added to the clique (having a None placeholder is helpful in indexing the weights of sigma)
        else:
            weight = 0
            for cnode in C['nodes']: 
                for edge in edges:
                    if node < cnode:
                        if [node, cnode] == edge[0:2]:
                            weight += edge[2] #Add the weight of each edge that connects the selected node to a node in the clique
                    elif cnode < node:
                        if [cnode, node] == edge[0:2]:
                            weight += edge[2] #Add the weight of each edge that connects the selected node to a node in the clique
            sigma.append(weight) 
    
    #Create I, T, and k
    I = [] #The mutually disjoint sets
    T = [] #The uncolored node set
    for snode in S:
        T.append(snode)
    k = 0 #The number of mutually disjoint sets
    
    #Create PI and upper
    PI = [] #The node sequence
    upper = [] #The upperbounds of each node in the sequence
    for weight in sigma:
        upper.append(0)
    
    #Continue until there are no more nodes left to color
    while T != []:
        
        #Create I[k] and X
        I.append([]) #Add Ik, which we initialize as an empty set
        X = [] #Candidate node set to add to I[k]
        for node in T:
            X.append(node)
        
        #Continue until there are no more nodes left to add to I[k]
        while X != []:
            
            #Find the node, v, in X with the smallest weight
            min_node = X[0]
            min_weight = sigma[X[0]]
            for node in X:
                if sigma[node] < min_weight:
                    min_node = node
                    min_weight = sigma[node]
            v = min_node
            
            I[k].append(v) #Add this v to I[k]
            
            PI = [v] + PI #Add v to the head of PI 
            
            #Update upper[v] by adding the sum of the max weights of a node given that that node is in a disjoint set smaller than the current disjoint set
            upper[v] += sigma[v]
            for i in range(k):
                max_weight = 0
                for node in I[i]:
                    if sigma[node] > max_weight:
                        max_weight = sigma[node]
                upper[v] += max_weight     
            
            X = [node for node in X if node not in N(v, partition_nodes)] #Remove nodes from X that are neighbors to v
            X.remove(v) #Remove v from X

            T.remove(v) #remove v from T
        
        #Look through every node that needs to be colored
        for v in T:
            
            #Update sigma by adding the max weight between v and a point in the disjoint set I[k] that is also a neighbor of v
            max_weight = 0
            for edge in edges:
                if v == edge[0]:
                    if edge[1] in I[k] and edge[1] in N(v, partition_nodes): #(note: in the case of a multipartite graphs, there is no need to check the neighbors of v)
                        if max_weight < edge[2]:
                            max_weight = edge[2] 
                if v == edge[1]:
                    if edge[0] in I[k] and edge[0] in N(v, partition_nodes): #(note: in the case of a multipartite graphs, there is no need to check the neighbors of v)
                        if max_weight < edge[2]:
                            max_weight = edge[2]
            sigma[v] += max_weight 
            
        k += 1 #Focus on next disjoint set
                 
    return PI, upper
                
#Find neighbors of a node using partition_nodes (note: this only works for complete multipartite graphs)
def N(v, partition_nodes):
    
    Neighbors = []
    for partition in partition_nodes:
        if v not in range(partition[0], partition[1]):
            Neighbors.extend(range(partition[0], partition[1])) #include all nodes in partition if partition doesn't contain the selected node
            
    return Neighbors

#Find neighbors of a node using edges        
def N_thorough(v, edges):
    
    Neighbors = []
    for edge in edges:
        if v == edge[0]:
            Neighbors.append(edge[1]) #Include all nodes that share an edge with the selected node
        elif v == edge[1]:
            Neighbors.append(edge[0]) #Include all nodes that share an edge with the selected node  
            
    return Neighbors
                
##############################################################################################################################################

#Return motif info from MECQ
def return_motif_info(cliques, all_motifs, seq_name_list):

    motif_solutions = []
    for Cmax in cliques:
        Cmax["nodes"].sort()
        motif_solution = {}
        seq_count = 0
        for node in Cmax["nodes"]: #Look through all solutions of one Cmax
            motif_solution[seq_name_list[seq_count]] = all_motifs[node] #Add selected motif to dict of one motif solution
            seq_count += 1
        score = Cmax["tot_weight"] #Return the score of the solution
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
                print("\n\nScore of solution = " + str(motif_solutions[n]["score"]) + "\n")
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
                print(motif_solutions[n]["score"], "\n")
            print("-" * 60)
            print("-" * 60)
            print("")

    #Print the motif solutions in a text file
    elif print_location == "txt_file":

        solution_filename = "B_and_B_method_Solutions_(" + filename.replace(".txt", "") + ")" + ".txt"
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
                print("\n\nScore of solution = " + str(motif_solutions[n]["score"]) + "\n", file = motif_file)
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
                print(motif_solutions[n]["score"], "\n", file = motif_file)
            print("-" * 60, file = motif_file)
            print("-" * 60, file = motif_file)
            print("\n\n\n\n\n", file = motif_file)

        motif_file.close()

    #Print the motif solutions in a csv file
    elif print_location == "csv_file":

        solution_filename = "B_and_B_method_Solutions_(" + filename.replace(".txt", "") + ")" + ".csv"
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
            score_row.append(motif_solution["score"])
        motif_writer.writerow(score_row)

        motif_writer.writerow([])
        motif_writer.writerow([])

        motif_file.close()
