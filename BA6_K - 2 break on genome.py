#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 01:49:48 2020

@author: jasonmoggridge

We can extend this pseudocode to a 2-break defined on genome P.

2-BreakOnGenome(P, i1 , i2 , i3 , i4 )
     GenomeGraph ← BlackEdges(P) and ColoredEdges(P)
     GenomeGraph ← 2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 )
     P ← GraphToGenome(GenomeGraph)
     return P

Code Challenge: Implement 2-BreakOnGenome.
Input: A genome P, followed by indices i1 , i2 , i3 , and i4 .
Output: The genome P' resulting from applying the 2-break operation
 2-BreakOnGenome(GenomeGraph i1 , i2 , i3 , i4 ).



Sample Input:

(+1 -2 -4 +3)
1, 6, 3, 8

Sample Output:

(+1 -2)(-3 +4)

"""
import copy
    


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
    
    print("\n\nChromosome to cycle:\n", Chromosome, '\nNodes ->\n', Nodes)
    return Nodes
#
    
def Cycle_to_Chromosome(Nodes):    
    
    Chromosome = [False for _ in range(0, len(Nodes),2)]
    for i in range( 0, len(Nodes),2):
        if Nodes[i] < Nodes[i+1]:
            Chromosome[i//2] = (Nodes[i]+1)//2
        else:
            Chromosome[i//2] = -int((Nodes[i]+1)//2)
    print('cycle to chromosome: \nNodes\n', Nodes,'\nChromosome:\n', Chromosome)
    return Chromosome
#


def Black_edges(genome):

    nodes= []
    for chromosome in genome:
        nodes = nodes + Chromosome_to_Cycle(chromosome)
    black_edges = []
    while nodes:
        black_edges.append((nodes.pop(0), nodes.pop(0), 'black'))
    print('\n\nBlack edges:\n', black_edges )
    return black_edges
#

def Colored_Edges(Chromosomes, colour):        

    col_edges = []
    for chromosome in Chromosomes:
        Nodes = Chromosome_to_Cycle(chromosome)
        for i in range(1, len(Nodes)-1, 2):
            col_edges.append((Nodes[i], Nodes[i+1], colour))
        col_edges.append((Nodes[-1],Nodes[0], colour))
    print('\n\nColoured edges:\n', col_edges)
    return col_edges
#

def Breakpoint_graph(edges):
    
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
#
    

def Two_Break_Genome(edges, i, colour):
    
    if (i[0], i[1], colour) in edges:
        edges.remove((i[0], i[1], colour))

    else:
        edges.remove((i[1], i[0], colour))
    if (i[2], i[3], colour) in edges:
        edges.remove((i[2], i[3], colour))
    else:
        edges.remove((i[3], i[2], colour))
    
    edges.append((i[0], i[2], colour))
    edges.append((i[1], i[3], colour))
    print('\n\nAfter breaking at ijkl', indices, '\nEdges:\n:', edges)
    return edges

#

# Dfs search function to enumerate cycles in breakpoint graph of 2 genomes

def DFS_cycles(breakpoint_graph):

    def dfs_explorer(v, colour, path, explored, cycle, graph): #
        
        print('\n\t dfs explorer: v', v, 'path:\n', path)
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
    
    print('\n\nDfs wrapper ---')
    graph = breakpoint_graph
    explored = dict(zip(graph.keys(), [False for _ in breakpoint_graph.keys()])) 
    cycles = [] 
    explored[1] = True
    dfs_explorer(2,'black',[1], explored, False, breakpoint_graph)
    
    #rest of chromosomes:
    for node in breakpoint_graph.keys():
        if not explored[node]: 
            dfs_explorer(node, 'red', [], explored, False, breakpoint_graph)     
    return cycles


### Main
   
#Parse -> genome -> 2-break (i,j,k,l)
    
with open("/Users/jasonmoggridge/Desktop/rosalind_ba6k.txt",'r') as infile:
    line = infile.readline().strip()
    line = line.strip()[1:-1]
    line = line.split(')(')
    genome = []
    for chromosome in line:
        genome.append([int(x) for x in chromosome.split(' ')])

    indices = infile.readline().strip()
    indices = [int(i) for i in indices.split(', ')]
    del (line, chromosome)
    infile.close()

# black edges represent the synteny blocks.
black_edges = Black_edges(genome)

#red edges represent the edges between synteny blocks
red_edges = Colored_Edges(genome,'red')

# Remove two edges and add two edges, specified by the two-break indices
edited_edges = Two_Break_Genome(red_edges, indices, 'red')

# Create an Adj list graph representation for rebuilding cycles
breakpointgraph = Breakpoint_graph(edited_edges + black_edges)

rebuilt_genome = DFS_cycles(breakpointgraph)




genome = []
for cycle in rebuilt_genome:
    chromosome = Cycle_to_Chromosome(cycle)
    string = '('
    for block in chromosome:
        if block > 0:
            string += '+' + str(block) + ' '
        else:
            string += str(block) + ' '
    string = string[:-1]+ ')'
    genome.append(string)

print(' '.join(genome))

