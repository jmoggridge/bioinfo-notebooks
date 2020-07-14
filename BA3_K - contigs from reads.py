#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 13:25:43 2019
@author: jasonmoggridge
http://rosalind.info/problems/ba3k/

BA3K:
    Generate Contigs from a Collection of Reads
    
    Maximum non-branching paths problem:
        
        A node v in a directed graph Graph is called a 1-in-1-out node if 
        its indegree and outdegree are both equal to 1, 
        i.e., in(v) = out(v) = 1.  
        We can rephrase the definition of a "maximal non-branching path"
        from the main text as a path whose internal nodes are 1-in-1-out
        nodes and whose initial and final nodes are not 1-in-1-out nodes.
        Also, note that the definition from the main text does not handle
        the special case when Graph has a connected component that is an
        isolated cycle, in which all nodes are 1-in-1-out nodes.

The MaximalNonBranchingPaths pseudocode below generates all non-branching
 paths in a graph. It iterates through all nodes of the graph that are not
 1-in-1-out nodes and generates all non-branching paths starting at each
 such node. 
 
 ***In a final step, MaximalNonBranchingPaths finds all isolated
 cycles in the graph.*** - I havent coded this in explicitly??? but correct?
 


    MaximalNonBranchingPaths(Graph)
        Paths ← empty list
        for each node v in Graph
            if v is not a 1-in-1-out node
                if out(v) > 0
                    for each outgoing edge (v, w) from v
                        Path ← the path consisting of single edge (v, w)
                        while w is a 1-in-1-out node
                            extend NonBranchingPath by the edge (w, u) 
                            w ← u
                        add NonBranchingPath to the set Paths
        for each isolated cycle Cycle in Graph
            add Cycle to Paths
        return Paths

"""



def DeBruijn_graph(kmers):        
    """Returns a DeBruijn graph for a set of overlapping kmers"""    
    graph = {}
    for kmer in kmers:
        if kmer[:-1] not in graph:
            graph[kmer[:-1]] = [kmer[1:]]
        else:
            graph[kmer[:-1]].append(kmer[1:])
    return graph
    #


    
def degree_array(graph):

    # set of all vertices 
    V = set(list(v for v in graph.keys()))
    for v in graph:
        V.update(graph[v])
    V = list(V)

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
        
    return deg
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
                    
                    path = str(v+w[-1])
                    while deg[w] == (1,1):
                        w = graph[w][0]
                        path += w[-1]
                    paths.append(path)
                    #add the path to the list of all before checking next edge from v
    return paths



###





f = open('//Users/jasonmoggridge/Desktop/dataset_205_5.txt', 'r')
reads = list(line.strip() for line in f.readlines())    

graph = DeBruijn_graph(reads)
deg = degree_array(graph)
paths = max_non_branching_paths(graph, deg)
    
with open('//Users/jasonmoggridge/Desktop/rosalind_ba3k_output.txt', 'w') as f:
    for path in paths:
        f.write(path+' ')
    f.close