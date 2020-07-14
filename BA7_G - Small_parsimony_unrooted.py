#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 08:37:05 2020

@author: jasonmoggridge
"""


# Simple Hamming distance function for edge scores
def HammingDistance(v, w):
    H=0
    for i in range(len(v)):
        if v[i] != w[i]:
            H +=1
    return H


### SMALL PARSIMONY FUNCTION FOR A SINGLE COLUMN IN ALIGNMENT
    
def SmallParsimony(tree, seq, n):
    
    # Scoring, pointers function (leafs -> root)    
    def SmallParsimony_Scoring(tree, column):
    
        """ Scoring and pointers from Tree: 
            strateies 
            -> dynamic programming > walk tree to root
            -> optimum parismony scores for each base, for each node, (for each column)
              memoize scores, pointers that yielded that score to retrace paths"""
                
    # Tag is 0 if scores for vertex need to be computed.

        root = 2*n-1
        Tag = [1 for _ in range(n)] + [0 for _ in range(root - n)]
#        print(root, Tag)

    # Memoization data structures:
        # init Scores dict S with -> s[leaf] = {bases:inf, seq[col]:0}
        # init Pointers dict P

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
#        print('initial S ->', S)
        for v in range(n,root):
            S[v] = {}
            
    # Keep computing scores until done entire tree -> condition: tag(v) = 1
    # update a 'ripe vertex': Tags(childs)=(1,1) but Tag(v) =0
        
        while 0 in Tag:
            
          # tree is numbered in reverese topological order so take first 0 in tag
            for v in range(n,root):
                if Tag[v] == 0:
                    son, daughter = tree[v][0], tree[v][1]
                    if Tag[son] and Tag[daughter]:
                        Tag[v] = 1
                        break
#            print('v',v)
          # Init scoring and pointer arrays for ripe vertex, find childs
            S[v] = {}
            P[v] = {}

#            print('son/daughter =', tree[v])
          # compute array of scores for vertex -> S[v] for each base k in alphabet
            for k in alpha:
                
             # init arrays for scores: S[son/daughter[base] + mis/match score
                son_scores = []
                daughter_scores = []
                for kk in alpha:
                    
                  # dist_a(k,kk) -> mismatch = 1, match = 0
                    a = 1
                    if k == kk:
                        a = 0
                  # arrays of scores for each base in S[child]
                    son_scores.append(S[son][kk] + a)
                    daughter_scores.append(S[daughter][kk] + a)
        
             # Take min weight edges from each child to create scores[v]
                S[v][k] = min(son_scores) + min(daughter_scores)
         
             # Store pointers son and daughter
                son_p = alpha[son_scores.index(min(son_scores))]
                daughter_p = alpha[daughter_scores.index(min(daughter_scores))]
                P[v][k] = (son_p, daughter_p)
        
#        print(S,'\n\n', P)
        return S, P
        ###

    
    def SmallParsimony_Backtracking(tree, S, P, seq, i):    
    
        """ 1 - > Select MinParsimony
            
                find base k with minParsimonyScoring at root v
                add base k to seq[root] at position i
                ready to reconstruct all seqs only using backpointers in P
        
            2 - > Rebuild seqs root-> /leafs
                iterate down from root (#2*n-1) to last internal v (#n).
                    store bases[v] into seq[v]
                    set bases[son] and [daughter] to Pointer[v][bases[v]]   """
            
        # find best scoring base at root. 
        # put that base as last element in array-> [bases]
            # initiates backwalking array [bases]
            
        best = float('inf')
        root = 2*n-2
        bases = [False for _ in range(2*n-1)]
        for k in S[root].keys():  
            if S[root][k] < best:           
                best = S[root][k]
                bases[root] = k   
#        print(best, bases)
        # Visit all nodes down from root to all parents of leaves.
            # update the bases for son, daughter from Pointers[node][base]      
            # add the base for the current node to ancestor sequence
        stack = [root]
        while stack:
            
            v = stack.pop(0)
            if v >= n:
#                print(v, bases[v], tree[v])
                k = bases[v]
                seq[v] += k       
                [son, daughter] = tree[v] 
                bases[son] = P[v][k][0]
                bases[daughter] = P[v][k][1]
                stack += [son, daughter]
        return seq   
    
    # Small Parsimony wrapper 
    
    parsimony_score = 0
    for i in range(len(seq[0])):
        S,P = SmallParsimony_Scoring(tree, i)
        parsimony_score += min(S[2*n-2].values())   
#        print(parsimony_score)
        seq = SmallParsimony_Backtracking(tree, S, P, seq, i) 
    
    return seq, parsimony_score
    ##
   
    
    
### MAIN


# Parsing input data 
    # n -> number leafs 
    # adjlist as edges for internal nodes [n....2n-1] -> v->w

#with open("dataset_10335_10.txt", 'r') as infile:
with open("rosalind_ba7g.txt", 'r') as infile:
    
    alpha = 'ACGT'
    tree = {}
    seq = []
    n = int(infile.readline().strip())

    for line in infile.readlines():
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

    seq += ['' for _ in range(n,2*n-1)]
    infile.close()
    del(v,w,line)
    #

# new_edges
root = len(tree.keys())
edge = tree[0][0]
tree[0] = [root]
tree[edge].append(root)
tree[edge].remove(0)
tree[root] = [0, edge]

# remove edges that are leaf-> root direction
stack = [root]
while stack:
    v = stack.pop(0)
    stack += tree[v]
    for w in tree[v]:
        tree[w].remove(v)
    

# Main function call on the now-rooted tree, sequences, # of leafs
seq, parsimony_score = SmallParsimony(tree, seq, n)

# remove stupid root edges and replace with edge that has total weight of both

tree[root] = []
tree[edge].append(0)


# Enumerate all edges for output lines -> 'v->w:distance'
    # keep track of edge distances to make sure seqs correspond
    # to the most Parsimonious tree -> dists == parsimony_score, else error somewhere

edges = []
dists = 0
for v in tree.keys():
    for w in tree[v]:
        dist = HammingDistance(seq[v], seq[w])
        dists += dist
        edges.append(seq[w] + '->' + seq[v] + ':' + str(dist))
        edges.append(seq[v] + '->' + seq[w] + ':' + str(dist))



# Sanity check if backtracking created the proper tree:
print(dists, '==', parsimony_score, 'else something gone awry')

# Writing output 
with open("smallparsimonyout3.txt", 'w') as outfile:
    outfile.write(str(parsimony_score))
    for edge in edges:
        outfile.write('\n')
        outfile.write(edge)
    outfile.close()