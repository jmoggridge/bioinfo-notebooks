#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 05:23:34 2020
@author: jasonmoggridge

BA6L: 2-Break Sorting Problem

    Input: Two genomes with circular chromosomes
        on the same set of synteny blocks.
        
    Output: The sequence of genomes resulting from
        applying a shortest sequence of 2-breaks 
        transforming one genome into the other.
    
2-break algo pseudocode

ShortestRearrangementScenario(P, Q)
     output P
     RedEdges ← ColoredEdges(P)
     BlueEdges ← ColoredEdges(Q)
     BreakpointGraph ← the graph formed by RedEdges and BlueEdges
     while BreakpointGraph has a non-trivial cycle Cycle
          (i2,i3)<-An arbitrary edge from BlueEdges in a non trivial red-blue cycle
          (i1,i2)<-An edge from RedEdges originating at node i1
          (i3,i4)<-an edge from RedEdges originating at node i3
          RedEdges ← RedEdges with edges (i1, i2) and (i3, i4) removed
          RedEdges ← RedEdges with edges (i2, i3) and (i4, i1) added
          BreakpointGraph ← the graph formed by RedEdges and BlueEdges
          P ← 2-BreakOnGenome(P, i1 , i3 , i2 , i4 )
          output P



"""

import copy


# Turn a chromosome into list of 2n nodes 
# (direction implied by edge intergers diff =+/-1)
def Chromosome_to_Cycle(Chromosome):
    
    Nodes =[False for _ in range(len(Chromosome)*2)]
    for j in range(1,len(Chromosome)+1):
        if Chromosome[j-1] > 0:
            Nodes[2*j-2] = 2*Chromosome[j-1] - 1
            Nodes[2*j-1] = 2*Chromosome[j-1]
        else:
            Nodes[2*j-2] = -2*Chromosome[j-1] 
            Nodes[2*j-1] = -2*Chromosome[j-1] - 1
    return Nodes
#

def Cycle_to_Chromosome(Nodes):    
    
    Chromosome = [False for _ in range(0, len(Nodes),2)]
    for i in range( 0, len(Nodes),2):
        if Nodes[i] < Nodes[i+1]:
            Chromosome[i//2] = (Nodes[i]+1)//2
        else:
            Chromosome[i//2] = -int((Nodes[i]+1)//2)
    return Chromosome
#


# Creates edges for synteny blocks to return genome representation from Adj list
def getBlackEdges(genome):

    nodes= []
    for chromosome in genome:
        nodes = nodes + Chromosome_to_Cycle(chromosome)
    black_edges = []
    while nodes:
        black_edges.append((nodes.pop(0), nodes.pop(0), 'black'))
    return black_edges

# Creates list of edges joining the ends of synteny blocks, for a genome
def getColoredEdges(genome, colour):        

    col_edges = []
    for chromosome in genome:
        Nodes = Chromosome_to_Cycle(chromosome)
        for i in range(1, len(Nodes)-1, 2):
            col_edges.append((Nodes[i], Nodes[i+1], colour))
        col_edges.append((Nodes[-1],Nodes[0], colour))
    print('\n\nColoured edges:\n', col_edges)
    return col_edges


# Creates and adjacency list graph for an edgelist -> for dfs and returning rebuilt genome
def getAdjList(edges):
    
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


# Function to find all cycles in the breakpoint graph, walking alternate coloured edges
def getCycles(breakpoint_graph):

    #wrapper function - standard dfs -> with coloured edges to stop back travel to parent
    def dfs_cycles(v, colour, path, explored, cycle, graph): #
        
        explored[v] = True
        path += [v]    
        for dest in graph[v]:                  
            u = dest[0]; col = dest[1]      
            if col != colour:
                if u == path[0]:
                    cycles.append(copy.deepcopy(path))
                    explored[u] = True
                    return
                else:
                    if not explored[u]:
                        dfs_cycles(u, col, path, explored, cycle, graph)             
                        if cycle:
                            return
        path.remove(v)       

    ### dfs wrapper main    
        # don't start off at node 1 for this breakpoint graph
    explored = dict(zip(breakpoint_graph.keys(), [False for _ in breakpoint_graph.keys()])) 
    cycles = [] 
    for node in breakpoint_graph.keys():
        if not explored[node]: 
            dfs_cycles(node, 'blue', [], explored, False, breakpoint_graph)     
    print('dfs -blue/red -> all cycles - >\n\t', cycles)
    return cycles

   
# del two edges (i0,i1)(i2,i3) replace with (i0,i2)(i1,i3)
    
def doTwo_Break(edges, i, colour):
    
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
    print('\n\nAfter breaking at ijkl', i, '\nEdges:\n:', edges)
    return edges



def Graph_To_Genome(graph): # takes edgelist, returns genome format?

    cycles = []
    cycle = []
    while graph:
        edge = graph.pop(0)
        if not cycle:
            if edge[0] % 2 == 0:
                end = edge[0] - 1
            else:
                end = edge[0] + 1
        cycle += [e for e in edge]
        if cycle[-1] == end:
            cycles.append([cycle[-1]]+cycle[:-1])
            cycle = []

    Genome = []
    for cycle in cycles:
        Genome.append(Cycle_to_Chromosome(cycle))
    return Genome




def getGenomeP(graph):

    def dfs_explorer(v, colour, path, explored, cycle, graph): #
        
        print('\n\t dfs explorer: v', v, 'path:\n', path)
        print('colour:', colour)
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
                        print('looking at u ->', u)
                        dfs_explorer(u, col, path, explored, cycle, graph)             
                        if cycle:
                            return
        # pop node off path
        path.remove(v)       

    # dfs wrapper main    
    
    print('\n\nDfs wrapper ---')
    explored = dict(zip(graph.keys(), [False for _ in graph.keys()])) 
    cycles = [] 
    explored[1] = True
    dfs_explorer(2,'black',[1], explored, False, graph)
    
    #rest of chromosomes:
    for node in graph.keys():
        if not explored[node]: 
            dfs_explorer(node, 'red', [], explored, False, graph)     
    return cycles


def display_new_P(P):
    
    genome = []
    for cycle in P:
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


### Main

# Parse genomes to 2 lists of lists of synteny blocks [[chromosome]...]
 
with open("/Users/jasonmoggridge/Desktop/rosalind_test.txt",'r') as infile:
    genomes = []
    for line in infile.readlines():
        line = line.strip()[1:-1]
        line = line.split(')(')
        genome =[]
        for chromosome in line:
            genome.append([int(x) for x in chromosome.split(' ')])
        genomes.append(genome)

    P,Q = genomes[0], genomes[1]
    del (line, chromosome, genomes)
    infile.close()

#

"""
ShortestRearrangementScenario(P, Q)
     output P
     RedEdges ← ColoredEdges(P)
     BlueEdges ← ColoredEdges(Q)
     BreakpointGraph ← the graph formed by RedEdges and BlueEdges
     while BreakpointGraph has a non-trivial cycle Cycle
          (i2,i3)<-An arbitrary edge from BlueEdges in a non trivial red-blue cycle
          (i1,i2)<-An edge from RedEdges originating at node i1
          (i3,i4)<-an edge from RedEdges originating at node i3
          RedEdges ← RedEdges with edges (i1, i2) and (i3, i4) removed
          RedEdges ← RedEdges with edges (i2, i3) and (i4, i1) added
          BreakpointGraph ← the graph formed by RedEdges and BlueEdges
          P ← 2-BreakOnGenome(P, i1 , i3 , i2 , i4 )
          output P """


BlackEdges = getBlackEdges(P)   
RedEdges = getColoredEdges(P, 'red')
BlueEdges = getColoredEdges(Q, 'blue')
BreakpointGraph = getAdjList(RedEdges + BlueEdges)


# Setup for while P != Q

blocks = 0
for p in P:
    blocks += len(p)
    
cycles = getCycles(BreakpointGraph) 

#while len(cycles) < blocks:
for cycle in  cycles:
    print('cycle - >', cycle)
    if len(cycle)>2:
        indices = [cycle[0], cycle[1], cycle[3], cycle[2]]
        RedEdges = doTwo_Break(RedEdges, indices, 'red')    
        Breakpointgraph = getAdjList(RedEdges + BlueEdges)
        P = getGenomeP(getAdjList(BlackEdges+RedEdges))
        
        print('P ->', P)
        break
        
        
        
        
        
        
        
        
        