#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 26 21:23:36 2020

@author: jasonmoggridge

Created on Sun Jan 26 20:12:14 2020

@author: jasonmoggridge

NeighborJoining(D)
    
    n ← number of rows in D
    
    if n = 2
        T ← tree consisting of a single edge of length D1,2
        return T
    
    D* ← neighbor-joining matrix constructed from the distance matrix D
    
    find elements i and j such that D*i,j is a minimum non-diagonal element of D*
    Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
    limbLengthi ← (1/2)(Di,j + Δ)
    limbLengthj ← (1/2)(Di,j - Δ)
    
    add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
    D ← D with rows i and j removed
    D ← D with columns i and j removed
    
    T ← NeighborJoining(D)
    
    add two new limbs (connecting node m with leaves i and j) to the tree T
    assign length limbLengthi to Limb(i)
    assign length limbLengthj to Limb(j)
    
    return T
    
    
# Stack of remaining clusters -> disambiguates cluster numbers from their D indices
    # eg. cluster name (ie Tree[cluster]) for row i in D = stack[i]
    # need to kep track of names of remaining leafs in D else get mixed up in recursion
    
        # eg. after join (0,1) to (5) from nrange -> [0,1,2,3,4]
        # D indexes[0..(n-1)] and stack = [2,3,4,5]

"""





def NeighborJoining(D, stack, new):
    
    # Start of Loop -> find out n nodes in matrix
    
    
    # setup D*, the NJ transformed matrix for neighbor joining min dist clusters
    n = np.size(D,1)
    D_star = np.zeros((n,n))
    new += 1
    # Intro header output
    print('\n>>> START OF NJ for D:\n at n->',n)    
    print('stack:', stack)
    print('new:', new)
    
    # Base Case: if n = 2, single edge left in matrix (dont join, unrooted tree)
    #   return ← tree consisting of a single edge of length D1,2

    if n == 2:            
        Tree = {}
        Tree[stack[0]] = [(stack[1], D[0][1])]
        Tree[stack[1]] = [(stack[0], D[0][1])]
        print('\n\t ###### THE BASE CASE -> single edge ##\nTree:\n', Tree)
        return Tree
        
    # Get total_dist on leaf x to all other leafs in matrix
    tot_dist = D.sum(axis=0)
        
    # Calculate all D*(i,j) elements using formula:
        # D*(i,j) = (n-2)* D(i,j) - sum(D[i]) - sum(D[j])    
    for i in range(n-1):
        for j in range(i+1,n):
            D_star[i][j] = (n-2)*D[i][j] - tot_dist[i] - tot_dist[j]
            D_star[j][i] = D_star[i][j]
    
    # Return (i,j) neighbors to join -> from argmin(D*)
    (i,j) = np.unravel_index(np.argmin(D_star, axis=None), D_star.shape) 
    print('\tArgmin(D_start) ->', (i,j))
    (Ci, Cj) = (stack[i], stack[j])
    
    # Calculate leaf edges to new parent the NJ way
    
    delta = (tot_dist[i] - tot_dist[j]) / (n-2)
    limb_i = (D[i][j] + delta)/2
    limb_j = (D[i][j] - delta)/2
    print('\t limb_i', limb_i, '\t limb_j', limb_j)
    del(delta, D_star, tot_dist)
    
    # Create array of new distances for new cluster[of i,j]
        # D[new][k] = (D[i][k] + D[j][k] - D[i][j]) / 2
        
    D_new = np.zeros((n-1))
    col = 0
    for k in range(n):
        if k not in (i,j):        
            D_new[col] = (D[i][k] + D[j][k] - D[i][j]) / 2
            col += 1
    
    
    # remove the rows and columns of D corresponding to Ci and Cj
    D = np.delete(D, (i, j), 0)
    D = np.delete(D, (i, j), 1) 
    
    
    # add a row to D for Cnew (without final diagonal entry...)
    D = np.append(D, [D_new[:-1]], 0)
    D = np.append(D, [[d] for d in D_new], 1)
    
    stack = np.append(stack, max(stack)+1)
    for leaf in (j,i):
        stack = np.delete(stack, leaf)
        
 
    # RECURSION STEP   
    Tree = NeighborJoining(D, stack, new)
    
    
    m = stack[-1]
    print('\tadding cluster m=',m, 'to tree...')
    Tree[m] += [(Ci, "%.3f" %(limb_i)), (Cj, "%.3f" %(limb_j))]
    Tree[Ci] = [(m, "%.3f" % (limb_i))]
    Tree[Cj] = [(m, "%.3f" % (limb_j))]
    print ('\n\n\tEnd of NJ function @n', n)
    return Tree




# Main #
 
import numpy as np

# Parse nodes n, n*n distance matrix D
with open("/Users/jasonmoggridge/Desktop/rosalind_test2.txt",'r') as infile:
    n = int(infile.readline())
    lines = [line.strip().split() for line in infile.readlines()]
    D = np.array(lines, dtype= float)
    infile.close()
    del(lines)
    
    new = n-1
    stack = np.array((range(n)))

Tree = NeighborJoining(D, stack, new) 

edges = []
for v in Tree.keys():
    for w in Tree[v]:
        edges.append((v, w[0], w[1]))

print('\n\n\n>>>NJ Tree for Distance Matrix D, as edgelist:\n\n')
# write edgelist to file
with open("/Users/jasonmoggridge/Desktop/NJ_tree_out.txt",'w') as outfile:       
    for edge in sorted(edges):
        print('\t{0}->{1}:{2}'.format(edge[0],edge[1],edge[2]))
        outfile.write('{0}->{1}:{2}\n'.format(edge[0],edge[1],edge[2]))
    outfile.close()