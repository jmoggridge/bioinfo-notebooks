#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 20:38:56 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3m/

BA3_M:
    Maximal Non-Branching Path Problem
        Find all maximal non-branching paths in a graph.
        
        Given: The adjacency list of a graph whose nodes are integers.
        
        Return: The collection of all maximal non-branching paths in the graph.

"""


def to_find(paths):
    found =set()
    for path in paths:
        found.update(path)
    
    return  V - found
#

    
def degree_array(graph):

    # set of all vertices 
    V = set(list(v for v in graph.keys()))
    for v in graph:
        V.update(graph[v])

    # compute deg for each v -> (deg_in, deg_out) 
    # all if/else statements -> looking up dict keys if/not created will throw except error

    deg = {}    
    for v in V:
            
        if v not in deg.keys():
            if v in graph.keys():
                deg[v] = (0, len(graph[v]))
            else:
                deg[v] = (0,0)
        else:
            if v in graph.keys():
                deg[v] = (deg[v][0], deg[v][1]+len(graph[v]))
        
        if v in graph.keys():
            for w in graph[v]:
                if w not in deg.keys():
                    deg[w] = (1,0)
                else:
                    deg[w] = (deg[w][0]+1,deg[w][1])
        
    return deg, V
    ##


    
def max_non_branching_paths(graph, deg):

    # list of all paths walked
    paths = []
    
    # need to check paths from all vertices
    for v in deg:
        
        # start walk if not at 1-in-1-out node
        if deg[v] != (1,1):
            
            # start a walk for all outgoing edges from branch node v
            if deg[v][1] > 0:
                
                for w in graph[v]:
                    # log the walk for each outgoing edge, start with 'v[:]+w[-1]'
                    # add last char for subseq nodes until reach another branch node
                    
                    path = [v,w]
                    while deg[w] == (1,1):
                        w = graph[w][0]
                        path.append(w)
                    paths.append(path)
                    #add the path to the list of all before checking next edge from v
    return paths
    #
    
def non_branching_cycles(graph):

    cycles=[]
    for v in not_found:
        if v not in found:
            found.update([v])
            if deg[v][1] > 0:
                for w in graph[v]:
                    cycle = [v]
                    while deg[w] == (1,1):
                        cycle.append(w)
                        w = graph[w][0]
                        if w == cycle[0]:
                            cycle.append(w)
                            cycles.append(cycle)
                            found.update(cycle)
                            break
    return cycles
    ##

##

#def non_branching_cycles

f = open('//Users/jasonmoggridge/Desktop/dataset_6207_2.txt', 'r')
lines = list(str(l.strip('\n')) for l in f.readlines())
Edges = [[i.strip(' ') for i in line.split('->')] for line in lines]

graph = {}
for edge in Edges:
    graph[int(edge[0])-1]=[int(i)-1 for i in edge[1].split(',')]
del((Edges, edge, lines, f))

deg, V = degree_array(graph)
paths = max_non_branching_paths(graph, deg)
not_found = 
cycles = non_branching_cycles(graph)



### Formatting output
adjusted =[]
for path in paths + cycles:
    adjusted.append(list(int(i+1) for i in path))

strings =[]    
for path in adjusted:
    strings.append( ' -> '.join(str(i) for i in path))
    
for string in strings:
    print(string)
    #print('\n\n\n',string)

#    ###

with open('//Users/jasonmoggridge/Desktop/rosalind_ba3m_output.txt', 'w') as f:
    for string in strings:
        f.write(string+'\n')
    f.close