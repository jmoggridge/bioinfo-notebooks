#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 00:49:11 2020

@author: jasonmoggridge

NearestNeighborInterchange(Strings)
    score ← ∞
    ->Small Parsimony in an Unrooted Tree Problem for Tree
        getTree(seqs)
        newScore ← the parsimony score of Tree
        newTree ← Tree
        newLabels ← labels on nodes of Tree
    
    while newScore < score
        score ← newScore
        Tree ← newTree
        for each internal edge e in Tree
            for each nearest neighbor NeighborTree of Tree with respect to the edge e
                solve the Small Parsimony in an Unrooted Tree Problem for NeighborTree
                neighborScore ← the minimum parsimony score of NeighborTree
                if neighborScore < newScore
                        newScore ← neighborScore
                        newTree ← NeighborTree
                        newLabels ← NeighborTreeLabels
        if newScore < score
            print newTree with newLabels

    return newTree


"""

import copy

# SUBROUTINES


def getTree(lines):

    seq = []
    tree = {}
    
    for line in lines:
        [v,w] = line.strip().split('->')
        if v[0] in alpha:
            if v in seq:
                v = seq.index(v)
            else:
                seq.append(v)
                v = len(seq)-1
        if w[0] in alpha:
            if w in seq:
                w = seq.index(w)
            else:
                seq.append(w)
                w = len(seq)-1
        v,w = int(v), int(w)
        if v not in tree.keys():
            tree[v] = [w]
        else:
            tree[v].append(w)
        if w not in tree.keys():
            tree[w] = [] 
    return tree, seq


# Root the tree:
def rootTree(tree):          
    
    rooted_tree = copy.deepcopy(tree)
    root = len(tree.keys())
    edge = tree[n][0]

    rooted_tree[n].append(root)
    rooted_tree[n].remove(edge)
    
    rooted_tree[edge].append(root)
    rooted_tree[edge].remove(n)
    
    rooted_tree[root] = [n, edge]
    
    stack = [root]
    while stack:
        v = stack.pop(0)
        stack += rooted_tree[v]
        for w in rooted_tree[v]:
            rooted_tree[w].remove(v)
    return rooted_tree


#   # Small Parsimony Scoring /Backtracking functions

def SP_Score(tree, column):
    
    root = 2*n-1
    Tag = [1 for _ in range(n)] + [0 for _ in range(root - n)]

    S = {}
    P = {} 
    for v in range(n):
        S[v] = {}
        P[v] = seq[v][column] 
        for k in alpha:
            if k == seq[v][column]:
                S[v][k] = 0
            else:
                S[v][k] = float('inf')

    for v in range(n,root):
        S[v] = {}
        
    while 0 in Tag:
        for v in range(n,root):
            if Tag[v] == 0:
                son, daughter = tree[v][0], tree[v][1]
                if Tag[son] and Tag[daughter]:
                    Tag[v] = 1
                    break
        S[v] = {}
        P[v] = {}
        for k in alpha:
            son_scores = []
            daughter_scores = []
            for kk in alpha:
                a = 1
                if k == kk:
                    a = 0
                son_scores.append(S[son][kk] + a)
                daughter_scores.append(S[daughter][kk] + a)
    
            S[v][k] = min(son_scores) + min(daughter_scores)
            son_p = alpha[son_scores.index(min(son_scores))]
            daughter_p = alpha[daughter_scores.index(min(daughter_scores))]
            P[v][k] = (son_p, daughter_p)

    return S, P


def SP_Seq(tree, S, P, seq, i):    

    best = float('inf')
    root = 2*n-2
    bases = [False for _ in range(2*n-1)]
    for k in S[root].keys():  
        if S[root][k] < best:           
            best = S[root][k]
            bases[root] = k   
    stack = [root]
    while stack:
        v = stack.pop(0)
        if v >= n:
            k = bases[v]
            seq[v] += k       
            [son, daughter] = tree[v] 
            bases[son] = P[v][k][0]
            bases[daughter] = P[v][k][1]
            stack += [son, daughter]
    return seq   

 
#  # Small Parsimony wrapper 
    
def SP_wrapper(tree, seq, n):
    
    rooted_tree = rootTree(tree)
    parsimony_score = 0
    seq = seq[:n] + ['' for _ in range(n,2*n-1)]

    for i in range(len(seq[0])):
        S,P = SP_Score(rooted_tree, i)
        parsimony_score += min(S[2*n-2].values())   
        seq = SP_Seq(rooted_tree, S, P, seq, i) 
   
    return seq, parsimony_score
####

def HammingDistance(v, w):
    
    H=0
    for i in range(len(v)):
        if v[i] != w[i]:
            H +=1
    return H

##

   

def NearestNeighbors(tree, edge):
    
    """  # use j,k indices to swap subtrees edges
      # creates two new trees that are -> nearest neighbors of tree, for edge (a,b)
      # b->j(z or y) and a->x. (a->w never breaks)           """
      
      
    # input internal edge to be broken for NN
    (a,b) = edge
    NN = [tree]
    
    # for wx|yz (tree0) Nearest Neighbors (x2) -> wy|xz (tree1), wz|xy (tree2) 
    for neigh in range(1,3):
        
        # copy original tree
        NN.append(copy.deepcopy(tree))
        
        # remove the a-b edge to easily get subtrees(a)->w,x and subtrees(b)-> yz
        NN[neigh][a].remove(b); NN[neigh][b].remove(a)
        [w,x] = NN[neigh][a]
        [y,z] = NN[neigh][b]

      # changes NN combo for the ->for NN loop to get the right change for each tree)
        if NN == 1:
            j, k = y, z
        else:
            j, k = z, y
   
     # break edge(b,j) and create edge(a,j)
      # edge(a,b) also added back to trees   
        NN[neigh][j].remove(b)
        NN[neigh][a] = [w,j,b]
        NN[neigh][j].append(a)   

      # break edge (a,x) and create edge(b,x)
        NN[neigh][x].remove(a)
        NN[neigh][b] = [x,k,a]
        NN[neigh][x].append(b)
        
    return NN[1], NN[2]
    ##


def writeTree(score, tree, seq, outfile):    
    # Enumerate all edges for output lines -> 'v->w:distance'
        # keep track of edge distances to make sure seqs correspond
        # to the most Parsimonious tree -> dists == parsimony_score, else error somewhere
    
    edges = []
    dists= 0
    for v in tree.keys():
        for w in tree[v]:
            dist = HammingDistance(seq[v], seq[w])
            dists+=dist
            edges.append(seq[v] + '->' + seq[w] + ':' + str(dist))
    
    print('dists ->', dists//2, 'vs score ->', score)
    # Writing output 
    outfile.write(str(score))
    for edge in edges:
        outfile.write('\n')
        outfile.write(edge)
    outfile.write('\n\n')



    
# MAIN
    
    
with open("rosalind_test.txt", 'r') as infile:
#with open("dataset_10336_8.txt", 'r') as infile:
    
    alpha = 'ACGT'
    n = int(infile.readline().strip())
    tree, seq = getTree(infile.readlines())
    infile.close()
    #
with open("large_parsimony_outfile.txt", 'w') as outfile:
    outfile.write('')


# Init -> score, given unrooted tree, SP_seqs, SP_score
score = float('inf')
newSeq, newScore = SP_wrapper(tree, seq, n)
newTree = tree
# iterating until Parsimony stops improving (lower score = better)

while newScore < score:
    # set current score, tree
    score = newScore
    tree = newTree

    internal = []
    # for all internal edges
    for a in tree.keys():
        if len(tree[a]) >1 :
            for b in tree[a]:
                if len(tree[b]) >1:
                    if (a,b) not in internal:
                        internal.append((a,b))

    for edge in internal:

        Neighbor_trees = NearestNeighbors(tree, edge[::-1])
        for neigh_tree in Neighbor_trees:
            neigh_labels, neigh_score = SP_wrapper(neigh_tree, seq, n)
            print('neigh_score ->', neigh_score)
            if neigh_score < newScore:
                newScore = neigh_score
                newTree = neigh_tree
                newSeq = neigh_labels

    if newScore < score:
        print('new score', newScore, '\n Neigh tree:\n', neigh_tree, newSeq)
        # Call output print function
        with open("large_parsimony_outfile.txt", 'a') as outfile:
            writeTree(newScore, newTree, newSeq, outfile)



