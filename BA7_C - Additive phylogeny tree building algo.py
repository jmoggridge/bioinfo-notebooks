#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:30:52 2020
@author: jasonmoggridge



AdditivePhylogeny(D)
    n ← number of rows in D
    if n = 2
        return the tree consisting of a single edge of length D1,2
    limbLength ← Limb(D, n)
    for j ← 1 to n - 1
        Dj,n ← Dj,n - limbLength
        Dn,j ← Dj,n
    (i, k) ← two leaves such that Di,k = Di,n + Dn,k
    x ← Di,n
    D ← D﻿ with row n and column n removed
    T ← AdditivePhylogeny(D)
    v ← the (potentially new) node in T at distance x from i on the path between i and k
    add leaf n back to T by creating a limb (v, n) of length limbLength
    return T



Code Challenge: Implement AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.

Input: An integer n followed by a space-separated n x n distance matrix.
Output: A weighted adjacency list for the simple tree fitting this matrix.



"""
######################################################################

### Additive Phylogeny subroutines


def Limb(D, n):
    """Return Limb_length(n) for current distance matrix 
        -> take mini over all (i,k) combinations in D_submatrix(:n)
            Limb(n) = min{ (D(i,n) + D(k,n) - D(i,k)) // 2}        """

    limb_len = float('inf') 
    for i in range(n-1):
        for k in range(n)[i+1:]:   

            if (D[i][n] + D[k][n] - D[i][k]) /2 < limb_len:
                limb_len = int((D[i][n] + D[k][n] - D[i][k]) /2)       
    return limb_len
##


def two_leaves_ik(D, n):
    """ Find (i,k) s.t. Limb(parent(n)) = 0 """
    
    for i in range(n-1):
        for k in range(i+1, n):
            if  D[i][k] == D[i][n] + D[k][n]:
                return (i, k)
##

def getPath(tree, source, dest):
    """Dfs wrapper: finds path(i -> k) in current tree
       -> returns list [(node, path_dist)..]    """
    
    
    def dfs(tree, source, dest):
        """iterative dfs over stack of nodes reachable from source"""
        
        stack = [(source, [source])]
        visited = set()
        
        while stack:
            (vertex, path) = stack.pop(0)
            if vertex not in visited:
                if vertex == dest:
                    return path
                visited.add(vertex)
                for neighbor in tree[vertex]:
                    stack.append((neighbor[0], path+[neighbor[0]]))
        return False
    
    # Wrapper
    path = dfs(tree, source, dest)
    
    # zip (nodes, cumulative path distances)
    path_d = [0]
    for i in range(1, len(path)):
        for dest in tree[path[i]]:
            if dest[0] == path[i-1]:
                path_d.append(path_d[i-1] + dest[1])

    return list(zip(path, path_d))
##



            
def AddLeaf(tree, path, n, limb, x):
    """ Adds next leaf and/or its parent to the expanding tree"""

    
    def addLeafParent(path, x):
    # Returns the parent node's index and updated tree with parent
        
        # algo: walks path(i,k) until path distance > or = x
            # case '=': vertex at x -> give index, same tree
            # case '>': no vertex at x, break enclosing edge of dist
                # add two new edges to/from parent in path(i,k)
        
        left = path[0]
        for right in path[1:]:
                
        # case path_d < D[i][n] : not far enough from source i  on path i -> k
            
            if right[1] < x:
                left = right
            
        # case '=': found exact vertex already in tree, ret: (parent, tree)
            
            elif right[1] == x:
                return right[0], tree
       
        # case '>': found appropriate edge to break if dist(path(i->right)) > distance_new node
            elif right[1] > x:
            
            # step1: break edge left -> right in tree (2 undir edges in dict/ for tree edge)
            
                tree[left[0]].remove((right[0], right[1] - left[1]))
                tree[right[0]].remove((left[0], right[1] - left[1]))
    
            
            # step2: calculate distance from left/right nodes to parent(n)
            
                left_x = x - left[1]
                right_x = right[1] - x
                
            
            # step3: add new node to tree with edges to left and right nodes
            
                new = len(tree)
                tree[new] = [(left[0], left_x), (right[0], right_x)]
                tree[left[0]].append((new, left_x))
                tree[right[0]].append((new, right_x))

                return new, tree        

    # Main tree-construction wrapper:
    
        # to add leaf n -> find where parent joins tree ()
        # depending if parent is already in tree 
            # iff parent already in tree, 
        # then add leaf n at parent(n)    
    
    parent_n, tree = addLeafParent(path, x)

    # adding leaf edge to parent
    tree[parent_n].append((n, limb))
    tree[n].append((parent_n, limb))
    return tree        


######################################################################


### Main wrapper function for recursive leaf joining for tree

def AdditivePhylogeny(D,n):
    
    # Base case: only two leafs -> return single edge between two nodes, no other steps    
    if n == 1:    
        nodes = list(range(len(D)))
        tree = dict(zip(nodes, list([] for _ in nodes)))
        tree[0].append((1, D[0][1]))
        tree[1].append((0, D[0][1]))
        return tree
    
    # (3) Create the D_bald matrix, where D[j][n] is d(j -> parent(n))   
    
    limb = Limb(D, n)
    for j in range(n):
        D[j][n] = D[j][n] - limb
        D[n][j] = D[j][n]
    
    # (4) Get (i, k) for limb_parent => min[(Din + Dkn - Dkn)/2]
    # x = dist(parent(n) -> i)
    
    (i,k) = two_leaves_ik(D, n)
    x = D[i][n]
    print('\n\tlimb(n) =', limb)
    print('\t(i,k) for n:', (i,k))
    print('\tx for n:', x)
    
    # (5) Trim distance matrix for recursive fx call

    tree = AdditivePhylogeny(D, n-1)
    
    print('\n\n Building tree now...')
  
    # (6) add leaf n and its parents to the tree on recursion up.

    path = getPath(tree,i,k)
    print('\tgot path: ->', path)
    
    tree = AddLeaf(tree, path, n, limb, x)
    print('\tTree build up to leaf:', n)

    return tree
#   
   
    

######################################################################
  

### Parse inputs:
    
#with open("/Users/jasonmoggridge/Desktop/rosalind_ba7c.txt",'r') as infile:
with open("/Users/jasonmoggridge/Desktop/rosalind_test1.txt",'r') as infile:

    n = int(infile.readline())

    D = []
    lines = infile.readlines()
    for line in lines:
        D.append(list(int(i) for i in line.strip().split()))
    del(lines, line)
    infile.close()
#


tree = AdditivePhylogeny(D, len(D)-1) 
    
lines = []
for key in tree.keys():
    for edge in tree[key]:
        string = str(key) + '->' + str(edge[0]) + ':' + str(edge[1])
        lines.append(string)

with open("/Users/jasonmoggridge/Desktop/additive_phylo_tree.txt",'w') as outfile:
    for line in lines:
        print(line)
        outfile.write(line+'\n')
outfile.close()
    