#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 22:53:06 2020
@author: jasonmoggridge

BA7D

UPGMA CLUSTERING

upgma psuedocode:
    
UPGMA(D, n):
      
Clusters ← n single-element clusters labeled 1, ... , n
     
 construct a graph T with n isolated nodes labeled by single elements 1, ... , n
    for every node v in T 
        Age(v) ← 0
    
    while there is more than one cluster
        find the two closest clusters Ci and Cj
        merge Ci and Cj into a new cluster Cnew with |Ci| + |Cj| elements

        add a new node labeled by cluster Cnew to T
        connect node Cnew to Ci and Cj by directed edges
 
        Age(Cnew) ← DCi, Cj / 2
        remove the rows and columns of D corresponding to Ci and Cj

        remove Ci and Cj from Clusters

        add a row/column to D for Cnew by computing D(Cnew, C) for each C in Clusters
        add Cnew to Clusters
    root ← the node in T corresponding to the remaining cluster
    for each edge (v, w) in T
        length of (v, w) ← Age(v) - Age(w)
    return T

* one issue to avoid is that the nodes in the distance matrix keep changing

soln -> keep a stack eg. (1,3,4,7,11,12) so that array indices match node in tree
        -> remove nodes from stack when they are clustered until on
"""

    
    
def UPGMA(D,n):    
    
 # init -> Clusters, Tree, Age, stack
    
    Clusters = {node:[node] for node in range(n)}
    Tree = {node:[] for node in range(n)}
    Age = {node:0 for node in range(n)}
    stack = list(Clusters.keys())

 # do while more than one cluster
    while len(stack)>1:

    # new cluster, index k
        Cnew = len(Tree)
        
    # find closest clusters i,j and get their distance
        [i, j] = np.unravel_index(np.argmin(D, axis=None), D.shape)    
        Ci, Cj = stack[i], stack[j]

    # Age(Cnew) ← Distance 0.5*D[i,j]
        Age[Cnew] = D[i][j]/2
        
    # merge Ci and Cj into a C_new / add Cnew to Clusters
        Clusters[Cnew] = Clusters[Ci] + Clusters[Cj]
        
    #add a new node Cnew to T, w directed edges to Ci and Cj
        Tree[Cnew] = [Ci, Cj]
        
    # computing D(Cnew, C) for each C in Clusters ->average_pairwise_distances
        
        size_Ci, size_Cj = len(Clusters[Ci]), len(Clusters[Cj])
        D_new = []
        
    # want to iterate over all columns m in matrix execpt i and j
        for m in range(len(stack)):
            if m != i and m != j:
                dist = ((D[i][m]*size_Ci) + (D[j][m]*size_Cj) ) / (size_Ci + size_Cj)
                D_new.append(dist)
        
    # remove Ci and Cj from stack
        stack.remove(Ci); stack.remove(Cj)
        
    # remove the rows and columns of D corresponding to Ci and Cj
        D = np.delete(D, (i, j), 0)
        D = np.delete(D, (i, j), 1) 
        
    # add a row to D for Cnew (without final diagonal entry...)
        D = np.append(D, [D_new], 0)

    # add diagonal for new cluster
        D_new.append(float('inf')) 
        D = np.append(D, [[d] for d in D_new], 1)
        
        # add new cluster to stack
        stack.append(Cnew)
    

    # Create edgelist
    # root ← the node in T corresponding to the remaining cluster   
    # for each edge (v, w) has length ← Age(v) - Age(w)
    edges = []
    for v in Tree.keys():
        for w in Tree[v]:
            edges.append((v, w, "%.3f" % (abs(Age[v] - Age[w]))))
            edges.append((w, v, "%.3f" % (abs(Age[v] - Age[w]))))
    return edges            

 # - end of UPGMA - #


# Main #
 
import numpy as np

# Parse nodes n, n*n distance matrix D
with open("/Users/jasonmoggridge/Desktop/rosalind_ba7d.txt",'r') as infile:
    n = int(infile.readline())
    lines = [line.strip().split() for line in infile.readlines()]
    D = np.array(lines, dtype= float)
    np.fill_diagonal(D, np.inf)
    infile.close()
    del(lines)

# generate edgelist of ultrametric tree by UPGMA   
edges = UPGMA(D,n)

# write edgelist to file
with open("/Users/jasonmoggridge/Desktop/UPGMA_tree_out.txt",'w') as outfile:       
    for edge in sorted(edges):
        print('{0}->{1}:{2}'.format(edge[0],edge[1],edge[2]))
        outfile.write('{0}->{1}:{2}\n'.format(edge[0],edge[1],edge[2]))
    outfile.close()
