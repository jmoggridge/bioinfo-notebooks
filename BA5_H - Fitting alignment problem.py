#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 17:03:31 2020
@author: jasonmoggridge


Code Challenge: Solve the Fitting Alignment Problem.

Input: Two nucleotide strings v and w, where v has length
 at most 1000 and w has length at most 100.
 
Output: A highest-scoring fitting alignment between v and w.
 Use the simple scoring method in which matches count +1 and
 both the mismatch and indel penalties are 1.


key changes to global -> fitting: gives free ride through first+last column.
    ie. can start+ end anywhere in longer seq v for no penalty
    -> finds optimal subsequence of v for v',w global alignment
    -> max score from last column. End alignment at start col.
    
    Implementation:    
    - init score matrix with all 0's in 0th-column, free taxi through seq v
    - init score matrix, for j in first row: -j
    * doesn't offer free taxi in w, just in v when i=0 or len(v) *
    
    -unlike finding max in entire matrix like in local alignment,
    find best scoring alignment by looking only at the last column of 
    scoring matrix
    Start alignment builder from this (i,m) node.
    
    Trim start of alignment to only the aligned portion of v by killing align builder 
    as soon as j=0. Prevents starting gaps in v.
    
    Good job dude
    
"""


def Fitting_Alignment(v,w):

    # Scoring is +1 for match -1 for sub or indel

    S = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]  
    S[0][0] = 0
    Pointers = [[-1 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
    Pointers[0][0] = 0
    for i in range(1, len(v) +1):
        S[i][0] = 0 
        Pointers[i][0] = 0
    for i in range(1, len(w) +1):      
        S[0][i] = -i
        Pointers[0][i] = 1
    
    # Fill out S matrix with max score values
    for i in range(1, len(v) +1):
        for j in range(1, len(w) +1):
            match = -1
            if v[i-1] == w[j-1]:
                match = 1
            S[i][j] = max(                
                     S[i -1][j] - 1,
                     S[i][j-1] - 1,
                     S[i-1][j-1] + match)
            
            # pointer code: down 0, right 1, diag 2
            if S[i][j] == S[i-1][j] - 1:
                Pointers[i][j] = 0
            elif S[i][j] == S[i][j-1] - 1:
                Pointers[i][j] = 1
            elif S[i][j] == S[i-1][j-1] - 1:
                Pointers[i][j] = 2
            elif S[i][j] == S[i-1][j-1] + 1:
                Pointers[i][j] = 2

    # Finding best scoring fitting alignment from the last column in S matrix
    # column is len(w) + 1; this finds the endee position if v, 
    # where the alignment can be started
    
    score = -float('inf')
    end = 0
    for i in range(len(v)+1):
        if S[i][len(w)] > score:
            score = S[i][len(w)]
            end = i
            
    return Pointers, end, score
#    

def Build_alignment(Pointers, v, w, end): 
    # walks through pointers in rev topo order to build aligned sequences

    align1 = ''
    align2 = ''
    i = end
    j = len(w)
    
    while j > 0: # while still bases in shorter seq
        
        if Pointers[i][j] == 0:            
            align1 += v[i-1]        
            align2 += '-'           
            i -= 1                  
        
        elif Pointers[i][j] == 1:          
            align1 += '-'
            align2 += w[j-1]
            j -= 1
        
        elif Pointers[i][j] == 2:
            align1 += v[i-1]
            align2 += w[j-1]
            i -= 1; j -= 1

    return align1[::-1], align2[::-1]
###


#with open("/Users/jasonmoggridge/Desktop/rosalind_ba5h.txt",'r') as infile:
#    v = infile.readline().strip()
#    w = infile.readline().strip()
#    infile.close()
#    
    
v = 'ACGACCACAGATACCGCTATTCACTATATCGTT'
w = 'GATACACT'

Pointers, end, score = Fitting_Alignment(v, w)
align = Build_alignment(Pointers, v, w, end)
#
print(score,'\n',align[0],'\n', align[1])

with open("/Users/jasonmoggridge/Desktop/rosalind_ba5h_out.txt",'w') as outfile:
    outfile.write(str(score) + '\n')
    outfile.write(align[0] + '\n')
    outfile.write(align[1])
    outfile.close()
