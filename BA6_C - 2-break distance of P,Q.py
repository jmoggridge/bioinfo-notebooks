#!/usr/bin/env python3
"""
Created on Fri Jan 17 20:28:36 2020

@author: jasonmoggridge
# -*- coding: utf-8 -*-


Code Challenge: Solve the 2-Break Distance Problem.

Input: Genomes P and Q.
Output: The 2-break distance d(P, Q).
Extra Dataset

You may be wondering how the graph representation that we have been using for breakpoint graphs could be transformed into an adjacency list. After all, we haven’t even labeled the nodes of this graph! Check out Charging Station: From Genomes to the Breakpoint Graph to see how to transform a genome into a graph that is easier to work with in our implementations.

Sample Input:

(+1 +2 +3 +4 +5 +6)
(+1 -3 -6 -5)(+2 -4)
"""



import copy
import sys
sys.setrecursionlimit(10**8)


# Functions


# returns number of shared blocks in P,Q

def Blocks(genomes):    
        
    blocks = []
    for chromosome in genomes[0]:
        for block in chromosome:
            blocks.append(abs(block))
    
    blocks_intersection = []
    for chromosome in genomes[1]:
        for block in chromosome:
            if abs(block) in blocks:
                blocks_intersection.append(block)
    return len(blocks_intersection)


# Parse synteny blocks from chromsome to nodes.

def Chromosome_to_Cycle(Chromosome):
    # Turns chromosome into list of 2n nodes (direction implied by edge intergers diff =+/-1)
    Nodes =[False for _ in range(len(Chromosome)*2)]
    for j in range(1,len(Chromosome)+1):
        if Chromosome[j-1] > 0:
            Nodes[2*j-2] = 2*Chromosome[j-1] - 1
            Nodes[2*j-1] = 2*Chromosome[j-1]
        else:
            Nodes[2*j-2] = -2*Chromosome[j-1] 
            Nodes[2*j-1] = -2*Chromosome[j-1] - 1
    return Nodes


# Edgelist of coloured edges from a genome

def Colored_Edges(Chromosomes, colour):        

    col_edges = []
    for chromosome in Chromosomes:
        Nodes = Chromosome_to_Cycle(chromosome)
        for i in range(1, len(Nodes)-1, 2):
            col_edges.append((Nodes[i], Nodes[i+1], colour))
        col_edges.append((Nodes[-1],Nodes[0], colour))
    return col_edges


# Create undirected breakpoint graph ->as adjacency list from both sets of coloured 
# edges

def Breakpoint_graph(genomes):
    
    edges= []
    colours = ['blue','red']
    for genome in genomes:
        edges += Colored_Edges(genome, colours.pop())
    Adj = {}
    for edge in edges:
        if edge[0] not in Adj.keys():
            Adj[edge[0]] = [(edge[1], edge[2])]
        else:  
            Adj[edge[0]].append((edge[1], edge[2]))
        if edge[1] not in Adj.keys():
            Adj[edge[1]] = [(edge[0],edge[2])]
        else:  
            Adj[edge[1]].append((edge[0],edge[2]))
    return Adj


# Dfs search function to enumerate cycles in breakpoint graph of 2 genomes

def DFS_cycles(breakpoint_graph):

    def dfs_explorer(v, colour, path, explored, cycle, graph): #
        
        # add node to current path and explored as True
        explored[v] = True; path += [v]    
    
        # dfs for all outgoing edges for v
        for dest in graph[v]:                  
            u = dest[0]; col = dest[1]      
            
            # avoid backedge to parent of v
            if col != colour:
                
                # cycle
                if u == path[0]:
                    cycles.append(copy.deepcopy(path))
                    explored[u] = True
                    return
            
                # not cycle
                else:
                    if not explored[u]:        
                        dfs_explorer(u, col, path, explored, cycle, graph)             
                        if cycle:
                            return
        # pop node off path
        path.remove(v)       

    # dfs wrapper main    
    graph = breakpoint_graph
    explored = dict(zip(graph.keys(), [False for _ in breakpoint_graph.keys()])) 
    cycles = []   
    for node in breakpoint_graph.keys():
        if not explored[node]: 
            path = []
            cycle = False
            dfs_explorer(node, 'blue', path, explored, cycle, breakpoint_graph)     

    return cycles


def Distance(genomes):
    return Blocks(genomes) - len(DFS_cycles(Breakpoint_graph(genomes)))

###
    # Main ---
        
with open("/Users/jasonmoggridge/Desktop/rosalind_ba6c.txt",'r') as infile:
    lines = infile.readlines()
    genomes =[]
    for line in lines:
        line = line.strip()[1:-1]
        line = line.split(')(')
        chromosomes = []
        for chromosome in line:
            chromosomes.append([int(x) for x in chromosome.split(' ')])
        genomes.append(chromosomes)
    infile.close(); del (line, lines, chromosomes) 
#
#blocks_pq = Blocks(genomes)
#breakpoint_graph = Breakpoint_graph(genomes)
#cycles = DFS_cycles(breakpoint_graph)


dist = Distance(genomes)
    
