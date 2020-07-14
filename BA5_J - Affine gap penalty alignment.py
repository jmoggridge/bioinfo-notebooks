#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 00:31:42 2020

@author: jasonmoggridge

BA5_J: Affine Gap Penalty Alignment

Alignment with Affine Gap Penalties Problem: Construct a
highest-scoring global alignment between two strings (with affine gap penalties).

Input:
    Two strings, a matrix score, and numbers σ and ε.
Output:
    A highest scoring global alignment between these strings,
    as defined by the scoring matrix score and by the gap opening
    and extension penalties σ and ε.


Keys to alog

-3 layer manhattan grid to store scores
- each grid has internal edges in only one direct + edges to adj grids up/down

- affine gap penalty p = sigma + k*epsilon
    where sigma is open penalty, epsilon is extend penalty, k is length of gap
- global alignment so start at 0,0 and end at n,m
- pointers matrix for each grid stores parent grid (indicates move, tf next position i,j)


"""


def Affine_Gap_Align(v,w):
    
    
    def Affine_gap_alignment_manhattan(v,w):
        
        #Global w affine gap penalty
    
    
        def build_Blosum62():
            """Build blosum62 mis/match scoring penalties mu"""
            
            with open("/Users/jasonmoggridge/Dropbox/Rosalind/textbook_track/Course3_Alignments/BLOSUM62.txt", 'r') as file:
                AA = file.readline().strip().split('  ')
                arrays = []
                for _ in AA:
                    array = list(file.readline().strip().split(' '))
                    while '' in array:
                        array.remove('')
                    array = list(int(i) for i in array[1:])
                    arrays.append(array)
                del(array)
    
            blosum62 = {}
            for i in range(len(AA)):
                for j in range(len(AA)):
                    blosum62[(AA[i],AA[j])] = arrays[i][j]
    
            return blosum62
    
        # Use the BLOSUM62 scoring matrix, a gap opening penalty of 11,
        # and a gap extension penalty of 1.
        blosum62 = build_Blosum62()
        sigma = 11
        epsilon = 1
    
        # Need three scoring matrices {low, mid, top} and three pointer matrices
        Low = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        Mid = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        Top = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        P_Low = [[False for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        P_Mid = [[False for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        P_Top = [[False for _ in range(len(w)+1)]for _ in range(len(v)+1)]
        
        # Initial row/col of 3 grids same scoring
        for grid in (Low, Mid, Top):
            for i in range(1, len(v) +1):
                grid[i][0] = -1* (sigma + (i-1) *epsilon)
            for j in range(1, len(w) +1):
                grid[0][j] = -1* (sigma + (j-1) *epsilon)
    
        # Init pointers grids, add backpointers to first row/col (low: up to i-1 row, top: j-1 col)
        for grid in (P_Low, P_Mid, P_Top):
            grid[0][0] = 0        
        for grid in (P_Low, P_Mid, P_Top):
            for i in range(1, len(v) +1):
                grid[i][0] = 'low'
            for j in range(1, len(w) +1):
                grid[0][j] = 'top'

        # Fill out 3 grids, Pointer P_[i][j] => grid of parent node, informs backedge move.
        
        for i in range(1, len(v)+1):
            for j in range(1, len(w)+1):
    
                # mu from blosum
                mu = blosum62[(v[i-1], w[j-1])]
                
                # max edge to low- either mid (open gap) or low (extend gap in w)
                Low[i][j] = max(Low[i-1][j] - epsilon,\
                          Mid[i-1][j] - sigma)
                if Low[i][j] == Low[i-1][j] - epsilon:
                    P_Low[i][j] = 'low'
                elif Low[i][j] == Mid[i-1][j] - sigma:
                    P_Low[i][j] = 'mid'
                
                # max edge to top - mid (open), top( ext. gap in v)
                Top[i][j] = max(Top[i][j-1] - epsilon,\
                          Mid[i][j-1] - sigma)
                if Top[i][j] == Top[i][j-1] - epsilon:
                    P_Top[i][j] = 'top'
                elif Top[i][j] == Mid[i][j-1] - sigma:
                    P_Top[i][j] = 'mid'
                
                # take max, take mis/match v[i],w[j] and move diag in mid 
                # or close gaps in either low(w) or top(v) and move up/down grids
                Mid[i][j] = max(Mid[i-1][j-1] + mu, Low[i][j], Top[i][j])
                if Mid[i][j] == Top[i][j]:
                    P_Mid[i][j] = 'top'
                elif Mid[i][j] == Low[i][j]:
                    P_Mid[i][j] = 'low'
                elif Mid[i][j] == Mid[i-1][j-1] + mu:
                    P_Mid[i][j] = 'mid'
        
        # need to take best score info to 
        best_score = max(Top[-1][-1], Mid[-1][-1], Low[-1][-1])
        return Low, Mid, Top, P_Low, P_Mid, P_Top, best_score
    ###
    
    def build_affine_gap_alignment(v,w,L,M,T,PL,PM,PT):
        
        V = ''
        W = ''
        i = len(v)
        j = len(w)
        
        for grid in (L,M,T):
            if grid[-1][-1] == score:
                start = grid
        if start == L:
            current = 'low'
        elif start == M:
            current = 'mid'
        elif start == T:
            current = 'top'
    
        while i >0 and j>0:
    
            if current == 'low':
                future = PL[i][j]
                V += v[i-1]
                W += '-'
                i -= 1
        
            if current == 'top':
                future = PT[i][j]
                V += '-'
                W += w[j-1]
                j -= 1
                    
            if current == 'mid':
                future = PM[i][j]
                if future == 'mid':
                    V += v[i-1]
                    W += w[j-1]
                    i -=1; j -= 1
    
            current = future
    
        return V[::-1], W[::-1]
    
    L, M, T, PL, PM, PT, score = Affine_gap_alignment_manhattan(v,w) 
    V,W = build_affine_gap_alignment(v,w,L,M,T,PL,PM,PT)
    return score,V, W
##
    

with open("/Users/jasonmoggridge/Desktop/rosalind_ba5j.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()
    infile.close()

score, V, W = Affine_Gap_Align(v,w)
print(score)
print(V)
print(W)