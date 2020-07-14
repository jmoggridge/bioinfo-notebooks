#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 00:48:34 2020

@author: jasonmoggridge


Given: Two protein strings s and t in FASTA format 
(with each string having length at most 1000 aa).

Return: The edit distance dE(s,t) followed by two 
augmented strings s′ and t′ representing an optimal alignment of s and t.
"""



def Edit_Alignment(v,w):
    S = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]  
    Pointers = [[-1 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
    S[0][0] = 0
    Pointers[0][0] = 0
    
    for i in range(1, len(v) +1):
        S[i][0] = -i         # need to penalize all the gaps in first col
        Pointers[i][0] = 0
    
    for i in range(1, len(w) +1):       # penalize all gaps in first row
        S[0][i] = -i
        Pointers[0][i] = 1
    
    
    for i in range(1, len(v) +1):
        for j in range(1, len(w) +1):
    
            match = -1
            if v[i-1] == w[j-1]:
                match = 0
    
            S[i][j] = max(                # pick longest path from 3 options
                     S[i -1][j] - 1,
                     S[i][j-1] - 1,
                     S[i-1][j-1] + match)
    
            if S[i][j] == S[i-1][j] - 1:
                Pointers[i][j] = 0
            elif S[i][j] == S[i][j-1] - 1:
                Pointers[i][j] = 1
            elif S[i][j] == S[i-1][j-1] - 1:
                Pointers[i][j] = 2
            elif S[i][j] == S[i-1][j-1]:
                Pointers[i][j] = 2
            
    return S, Pointers

def Build_alignment(Pointers, v, w): 
    # walks through pointers in rev topo order to build aligned sequences

    align1 = ''
    align2 = ''
    i = len(v)
    j = len(w)
    
    while i > 0 or j > 0: # different from LCS to finish global alignment
        
        if Pointers[i][j] == 0:            # pointer say go up a row
            align1 += v[i-1]        # add the base to align1  
            align2 += '-'           # add a gap to align2
            i -= 1                  # decrement row number
        
        elif Pointers[i][j] == 1:          # pointer says up move, same as above for seq b
            align1 += '-'
            align2 += w[j-1]
            j -= 1
        
        elif Pointers[i][j] == 2:          #diagonal move, respective bases to each align
            align1 += v[i-1]
            align2 += w[j-1]
            i -= 1; j -= 1
    
    #reverse path to have it start at first base for both seqs
    return align1[::-1], align2[::-1]

###

with open("/Users/jasonmoggridge/Desktop/rosalind_edta.txt",'r') as infile:
    
    infile.readline()
    lines = list(line.strip() for line in infile.readlines())
    seqs =[]
    seq = ''    
    for line in lines:
        if line[0] == '>':
            seqs.append(seq)
            seq = ''
        else:
            seq += str(line.strip())
    seqs.append(seq)
    
v, w = seqs[0], seqs[1]

S, Pointers= Edit_Alignment(v,w)
edit_distance = -S[-1][-1]
align = Build_alignment(Pointers, v,w)


with open("/Users/jasonmoggridge/Desktop/rosalind_edta_out.txt",'w') as outfile:
    outfile.write(str(edit_distance) + '\n')
    outfile.write(align[0] + '\n')
    outfile.write(align[1])
    outfile.close()
