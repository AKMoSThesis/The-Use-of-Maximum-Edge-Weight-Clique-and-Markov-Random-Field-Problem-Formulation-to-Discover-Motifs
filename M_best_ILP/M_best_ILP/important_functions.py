#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 12:59:07 2020

@author: iamserious
"""
import math
from gurobipy import *
from ortools.sat.python import cp_model
import csv
import matplotlib.pyplot as plt
import networkx as nx

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

######################################################### M_best_ILP-specific ############################################################

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

#Using Gurobi to solve ILP
def Gurobi_ILPmotif(num_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights):
    
    #Create model
    model = Model("poolsearch") 
    
    #Add decision variables x and y 
    x = model.addVars(node_count, vtype=GRB.BINARY, name="x") # Binary decision variable for the nodes
    y = model.addVars(edge_count, vtype=GRB.BINARY, name="y") # Binary decision variable for the edge between two nodes
    #print("set up variables DONE")
    
    #Set objective function
    model.setObjective(sum(weights[i] * y[i] for i in range(edge_count)), GRB.MAXIMIZE) 
    #print("set up objective function DONE")
    
    #Add node constraint 
    for j in range(partition_count):
        first_node, last_node = partition_nodes[j]
        model.addConstr(sum(x[u] for u in range(first_node, last_node)) == 1) 
    #print("set up node constraint DONE")
    #Add edge constraint 
    for v in range(node_count):
        model.addConstr(sum(y[i] for i in range(edge_count) if edges[i][0] == v or edges[i][1] == v) == (partition_count - 1) * x[v]) #(partition_count-1) because we are looking at all v in V\Vj
    #print("set up edge constraint DONE")
    
    #Specify desire for multiple solutions
    model.setParam(GRB.Param.PoolSolutions, num_solutions) #Specify number of desired solutions
    model.setParam(GRB.Param.PoolSearchMode, 2) #Specify desire to find n best solutions
    
#    #Write out the model
#    model.write("ILPmotif_model.lp")
    
    #Solve the model
    model.optimize()
    
    #Find selected nodes and edges of each solution
    ILPmotif_solutions = []
    for n in range(num_solutions):
        solution = []
        model.setParam(GRB.Param.SolutionNumber, n) #Look through each solution 
        for var in model.getVars(): #Look at each variable in one solution
            if var.Xn > 0.5: #Add selected node or edge of one solution
                solution.append(var.varName) 
        ILPmotif_solutions.append({"solution" : solution, "objective value" : model.PoolObjVal}) #Add solution and objective value to list of solutions
    
    return ILPmotif_solutions

#VarArraySolution taken from nb 3 solved 
class VarArraySolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, xvariables, yvariables, weights):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__xvariables = xvariables #ADDED
        self.__yvariables = yvariables #ADDED
        self.__weights = weights
        self.__solution_count = 0
        self.ILPmotif_solutions = [] #ADDED

    def OnSolutionCallback(self):
        self.__solution_count += 1
        solution = [] #ADDED
        var_count = 0 #ADDED
        for x in self.__xvariables:
#            print('%s=%i' % (v, self.Value(v)), end=' ') #REMOVED
            if self.Value(x) > 0.5: #ADDED
                solution.append('x['+str(var_count)+']') #ADDED
            var_count += 1 #ADDED
        var_count = 0 #ADDED
        objective_value = 0.0 #ADDED
        for y in self.__yvariables:
#            print('%s=%i' % (v, self.Value(v)), end=' ') #REMOVED
            if self.Value(y) > 0.5: #ADDED
                solution.append('y['+str(var_count)+']') #ADDED
                objective_value += self.__weights[var_count] #ADDED  
            var_count += 1 #ADDED
        self.ILPmotif_solutions.append({"solution" : solution, "objective value" : objective_value}) #ADDED
#        print() #REMOVE                    

    def SolutionCount(self):
        return self.__solution_count

#Using ORtools to solve ILP
def ORtools_ILPmotif(multiple_solutions, node_count, partition_count, edge_count, partition_nodes, edges, weights):
    
    #Finding the best solution
    if multiple_solutions == "False": 
        
        #Create model
        model1 = cp_model.CpModel()

        #Add decision variables x (node) and y (edge)
        x = [model1.NewIntVar(0, 1, 'x'+str(i)) for i in range(node_count)]
        y = [model1.NewIntVar(0, 1, 'y'+str(i)) for i in range(edge_count)]
        #print("set up variables DONE")

        #Add node constraint 
        for j in range(partition_count):
            first_node, last_node = partition_nodes[j]
            model1.Add(sum(x[u] for u in range(first_node, last_node)) == 1)
        #print("set up node constraint DONE")
        #Add edge constraint
        for v in range(node_count):
            model1.Add(sum(y[i] for i in range(edge_count) if edges[i][0] == v or edges[i][1] == v) == (partition_count - 1) * x[v])   
        #print("set up edge constraint DONE")
        
        #Set objective function (IF OPTIMAL SOLUTION DESIRED)
        model1.Maximize( sum(int(100000 * weights[i]) * y[i] for i in range(edge_count)) )
        #print("set up objective function DONE")

        #Solve the model
        solver1 = cp_model.CpSolver()
        status1 = solver1.Solve(model1)

        #Find selected nodes and edges of each solution
        ILPmotif_solutions = []
        solution = []
        for j in range(node_count): #Look through each node
            if solver1.Value(x[j]) > 0.5: #Add selected node of one solution
                solution.append('x['+str(j)+']')
        for j in range(edge_count): #Look through each edge
            if solver1.Value(y[j]) > 0.5: #Add selected edge of one solution
                solution.append('y['+str(j)+']')
        ILPmotif_solutions.append({"solution" : solution, "objective value" : solver1.ObjectiveValue()/100000}) #Add solution and objective value to list of solutions
        
        return ILPmotif_solutions
    
    #Finding multiple good solutions
    elif multiple_solutions == "True":
        
        #Find best solution
        ILPmotif_solutions = ORtools_ILPmotif("False", node_count, partition_count, edge_count, partition_nodes, edges, weights)
        best_objective_value = ILPmotif_solutions[0]["objective value"]
        
        #Create model
        model = cp_model.CpModel()

        #Add decision variables x (node) and y (edge)
        x = [model.NewIntVar(0, 1, 'x'+str(i)) for i in range(node_count)]
        y = [model.NewIntVar(0, 1, 'y'+str(i)) for i in range(edge_count)]
        #print("set up variables DONE")

        #Add node constraint 
        for j in range(partition_count):
            first_node, last_node = partition_nodes[j]
            model.Add( sum(x[u] for u in range(first_node, last_node)) == 1 )
        #print("set up node constraint DONE")
        #Add edge constraint
        for v in range(node_count):
            model.Add( sum(y[i] for i in range(edge_count) if edges[i][0] == v or edges[i][1] == v) == (partition_count - 1) * x[v] )   
        #print("set up edge constraint DONE")
        
        #ADD NEW CONSTRAINT
        new_const = model.Add( sum(int(100000 * weights[i]) * y[i] for i in range(edge_count)) >= int((0.99) * 100000 * best_objective_value))
        #print("set up new constraint DONE")
        
        #Solve the model
        solver = cp_model.CpSolver()
        solutionPrinter = VarArraySolutionPrinter(x, y, weights) #VarArraySolutionPrinterWithLimit is a class taken from nb 3 solved 
        status = solver.SearchForAllSolutions(model, solutionPrinter)        
        
        #Recover selected nodes and edges of each solution
        ILPmotif_solutions = solutionPrinter.ILPmotif_solutions
        
        return ILPmotif_solutions        

###########################################################################################################################################

#Return motif info from ILP solution
def return_motif_info(ILPmotif_solutions, all_motifs, seq_name_list): 
    
    motif_solutions = []
    for ILPmotif_solution in ILPmotif_solutions:
        motif_solution = {}
        seq_count = 0
        for var in ILPmotif_solution["solution"]: #Look through all solutions of ILPmotif
            if var[0] == "x": #Return the selected motif
                index = int(var.replace("x[", "").replace("]", ""))
                motif_solution[seq_name_list[seq_count]] = all_motifs[index] #Add selected motif to dict of one motif solution
                seq_count += 1 
        score = ILPmotif_solution["objective value"] #Return the score of the solution
        motif_solutions.append({"solution" : motif_solution, "score" : score})
        
    return motif_solutions

#Print the motif solutions in console
def print_motif_solutions(motif_solutions, solution_details, print_location = "console", print_format = "elegant"):
    
    #Unpack solution details
    filename, solver, motif_length, num_solutions, total_time = solution_details
    
    #Print the motif solutions in the console
    if print_location == "console":
        
        print("\nSOLUTIONS GENERATED!")
        print("filename =", filename)
        print("solver =", solver)
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
        
        solution_filename = "M_Best_ILP_Solutions_(" + solver + ")_(" + filename.replace(".txt", "") + ")" + ".txt"
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
        
        solution_filename = "M_best_ILP_Solutions_(" + solver + ")_(" + filename.replace(".txt", "") + ")" + ".csv"
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
     
#Return draw info from ILP solution
def return_draw_info(ILPmotif_solutions, edges):
    
    draw_solutions_nodes = []
    draw_solutions_edges = []
    for ILPmotif_solution in ILPmotif_solutions:
        draw_solution_nodes = []
        draw_solution_edges = []
        for var in ILPmotif_solution["solution"]: #Look through all solutions of ILPmotif
            if var[0] == "x": #Return the selected motif
                index = int(var.replace("x[", "").replace("]", ""))
                draw_solution_nodes.append(index) #Add selected node to list of motifs that will be drawn
            else:
                index = int(var.replace("y[", "").replace("]", ""))
                draw_solution_edges.append(edges[index])
        draw_solutions_nodes.append(draw_solution_nodes)
        draw_solutions_edges.append(draw_solution_edges)
        
    return draw_solutions_nodes, draw_solutions_edges
    
#Draw each solution as a max-weight clique and save as a separate png file          
def draw_graph(draw_solutions_nodes, draw_solutions_edges, partition_count, edge_count, partition_nodes, edges):  
    
    for n in range(len(draw_solutions_nodes)):
        
        #Select the edges and edges of one solution
        draw_solution_nodes = draw_solutions_nodes[n]
        draw_solution_edges = draw_solutions_edges[n]
    
        # Make graph
        G = nx.Graph()

        # Add nodes to graph
        for i in range(partition_count):
            first_node, last_node = partition_nodes[i]
            y_position = 0
            for node in range(first_node, last_node):
                y_position += 1
                for draw_node in draw_solution_nodes:
                    if draw_node == node:
                        G.add_node(node, pos=(i, y_position), node_color='pink') #Color selected node as pink
                        break
                    else:
                        G.add_node(node, pos=(i, y_position), node_color='skyblue') #Color non-selected node as blue

        # Add edges to graph     
        for i in range(edge_count):
            for j in range(len(draw_solution_edges)):
                if edges[i] == draw_solution_edges[j]:
                    G.add_edge(edges[i][0], edges[i][1], edge_color='pink') #Color selected edge as pink
                    break
                else:
                    G.add_edge(edges[i][0], edges[i][1], edge_color='none') #Color non-selected edge as none

        #Assemble dictionaries of all node and edge attributes
        node_position_dict = nx.get_node_attributes(G,'pos') 
        node_color_dict = nx.get_node_attributes(G,'node_color')
        edge_color_dict = nx.get_edge_attributes(G,'edge_color')
    
        #Assemble lists of node and edge colors
        node_color_list = node_color_dict.values()
        edge_color_list = edge_color_dict.values()

        #Draw graph
        nx.draw(G, pos = node_position_dict, node_color = node_color_list, edge_color = edge_color_list, with_labels = True) # Draw the graph
        plt.savefig("GRAPH_"+str(n)+".png")
        
        #Clear previous graph
        plt.clf()
    