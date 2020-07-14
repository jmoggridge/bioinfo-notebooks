#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 01:46:52 2020
@author: jasonmoggridge

BA7-F: SMALL PARSIMONY ALGORITHM

    - Return sequence labels for all ancestor nodes, Parsimony
    - Given tree and leaf seqs
    
    strategy
        - treat each column as independent subproblem
        - use dynamic programming to calc the recurrence relation
            - add one to score if mismatch between v->son/daughter -> 'a'
            Score[v][k] =   min { S[son][i] + a(i,k) } 
                            min{ S[daughter][j] + a(j,k) }
            Pointers[v][k] => (argmin_k(S[son] ), argmin_k(S[daughter])

PSEUDOCODE from rosalind BA7_f:
    
SmallParsimony(T, Character)
    for each node v in tree T
        Tag(v) ← 0
        if v is a leaf
            Tag(v) ← 1
            for each symbol k in the alphabet
                if Character(v) = k
                    sk(v) ← 0
                else
                    sk(v) ← ∞
    while there exist ripe nodes in T
        v ← a ripe node in T
        Tag(v) ← 1
        for each symbol k in the alphabet
            sk(v) ← minimumall symbols i {si(Daughter(v))+αi,k} + minimumall symbols j {sj(Son(v))+αj,k}
    return minimum over all symbols k {sk(v)}

*Doesn't say anything about how to walk backpointers ^ 
"""

# Simple Hamming distance function for edge scores
def HammingDistance(v, w):
    H=0
    for i in range(len(v)):
        if v[i] != w[i]:
            H +=1
    return H


### SMALL PARSIMONY FUNCTION FOR A SINGLE COLUMN IN ALIGNMENT
    
def SmallParsimony(Tree, seq, n):
    
    # Scoring, pointers function (leafs -> root)    
    def SmallParsimony_Scoring(Tree, column):
    
        """ Scoring and pointers from Tree: 
            strateies 
            -> dynamic programming > walk tree to root
            -> optimum parismony scores for each base, for each node, (for each column)
              memoize scores, pointers that yielded that score to retrace paths"""
                
    # Tag is 0 if scores for vertex need to be computed.

        root = 2*n-1
        Tag = [1 for _ in range(n)] + [0 for _ in range(root - n)]
        

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
        
        
    # Keep computing scores until done entire tree -> condition: tag(v) = 1
    # update a 'ripe vertex': Tags(childs)=(1,1) but Tag(v) =0
        
        while 0 in Tag:
            
          # tree is numbered in reverese topological order so take first 0 in tag
            for v in range(root):
                if Tag[v] == 0:
                    Tag[v] = 1
                    break
            
          # Init scoring and pointer arrays for ripe vertex, find childs
            S[v] = {}
            P[v] = {}
            son, daughter = Tree[v][0], Tree[v][1]
            
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
        
        return S, P
        ###

    
    def SmallParsimony_Backtracking(Tree, S, P, seq, i):    
    
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
        
        # Visit all nodes down from root to all parents of leaves.
            # update the bases for son, daughter from Pointers[node][base]      
            # add the base for the current node to ancestor sequence
    
        for v in range(2*n-2, n-1, -1):    
            k = bases[v]
            seq[v] += k       
            [son, daughter] = Tree[v] 
            bases[son] = P[v][k][0]
            bases[daughter] = P[v][k][1]

        return seq   
    
    # Small Parsimony wrapper 
    
    parsimony_score = 0
    for i in range(len(seq[0])):
        S,P = SmallParsimony_Scoring(Tree, i)
        parsimony_score += min(S[2*n-2].values())    
        seq = SmallParsimony_Backtracking(Tree, S, P, seq, i) 
    
    return seq, parsimony_score
    ##
    
    
    
### MAIN


# Parsing input data 
    # n -> number leafs 
    # adjlist as edges for internal nodes [n....2n-1] -> v->w

with open("dataset_10335_10.txt", 'r') as infile:
#with open("rosalind_test.txt", 'r') as infile:
    
    alpha = 'ACGT'
    Tree = {}
    seq = []
    

    n = int(infile.readline().strip())
    for line in infile.readlines():
        [v,w] = line.strip().split('->')
        
        if w[0] in alpha:
            seq.append(w)
            w = len(seq) -1
    
        v,w = int(v),int(w)
        
        if v not in Tree.keys():
            Tree[int(v)] = [w]
        else:
            Tree[v].append(w)
        if w not in Tree.keys():
            Tree[w] = [] 
    seq += ['' for _ in range(n,2*n-1)]
    infile.close()
    del(v,w, line)
    #
    
# Main function call on tree, sequences, # of leafs
seq, parsimony_score = SmallParsimony(Tree, seq, n)


# Enumerate all edges for output lines -> 'v->w:distance'
    # keep track of edge distances to make sure seqs correspond
    # to the most Parsimonious tree -> dists == parsimony_score, else error somewhere

edges = []
dists = 0
for v in Tree.keys():
    for w in Tree[v]:
        dist = HammingDistance(seq[v], seq[w])
        dists += dist
        edges.append(seq[w] + '->' + seq[v] + ':' + str(dist))
        edges.append(seq[v] + '->' + seq[w] + ':' + str(dist))


# Sanity check if backtracking created the proper tree:
print(dists, '==', parsimony_score, 'else something gone awry')

# Writing output 
with open("smallparsimonyout2.txt", 'w') as outfile:
    outfile.write(str(parsimony_score))
    for edge in edges:
        outfile.write('\n')
        outfile.write(edge)
    outfile.close()