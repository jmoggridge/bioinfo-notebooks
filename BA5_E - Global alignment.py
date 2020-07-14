#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 17:14:00 2020

@author: jasonmoggridge

Code Challenge: Solve the Global Alignment Problem.

Input: Two protein strings written in the single-letter
 amino acid alphabet.

Output: The maximum alignment score of these strings followed 
    by an alignment achieving this maximum score. Use the BLOSUM62
    scoring matrix for matches and mismatches as well as the indel
    penalty σ = 5.
"""

#import sys
#sys.setrecursionlimit(10**4)
#
def global_alignment_BLOSUM62(a,b):

    
    """ Create BLOSUM62 aa -> aa mis/match scoring matrix 
    # as a dictionary[aa,aa]= score"""
    
    def Blosum62():
        
        with open("/Users/jasonmoggridge/Desktop/BLOSUM62.txt", 'r') as file:
            AA = file.readline().strip().split('  ')
            arrays = []
            for _ in AA:
                array = list(file.readline().strip().split(' '))
                while '' in array:
                    array.remove('')
                array = list(int(i) for i in array[1:])
                arrays.append(array)
            del(array)
        
        BLOSUM62 = {}
        for i in range(len(AA)):
            for j in range(len(AA)):
                BLOSUM62[(AA[i],AA[j])] = arrays[i][j]
        return BLOSUM62    
    
    """ This function walks the manhattan grid and stores the best score to 
    # point (i,j) in the Scores matrix """
    
    def Scoring_and_Backtracks(a,b):
        
        Scores = [[0 for _ in range(len(b)+1)]for _ in range(len(a)+1)]  
        Backtracks = [[-1 for _ in range(len(b)+1)]for _ in range(len(a)+1)]
        
        Scores[0][0] = 0
        Backtracks[0][0] = 0
        
        for i in range(1, len(a) +1):
            Scores[i][0] = -i * sigma         # need to penalize all the gaps in first col
            Backtracks[i][0] = 0
            
        for i in range(1, len(b) +1):       # penalize all gaps in first row
            Scores[0][i] = -i * sigma
            Backtracks[0][i] = 1
        
        for i in range(1, len(a) +1):
            for j in range(1, len(b) +1):
    
                match = BLOSUM62[(a[i-1], b[j-1])]
                
                score = max(                # pick longest path from 3 options
                 Scores[i -1][j] - sigma,
                 Scores[i][j-1] - sigma,
                 Scores[i-1][j-1] + match)
                Scores[i][j] = score
                
                if score == Scores[i][j-1] - sigma:
                    Backtracks[i][j] = 1
                elif score == Scores[i-1][j] - sigma:
                    Backtracks[i][j] = 0
                elif score == Scores[i-1][j-1] + match:
                    Backtracks[i][j] = 2
                    
        return Scores, Backtracks
     
    
    """ Function rebuilds the alignment by walking along
    backpointers through the grid, adding a base or gap as indicated by 
    [down: 0, right: 1, diag: 2]    """
    
    def Build_alignment(Backtracks, a, b): 
        # walks through pointers in rev topo order to build aligned sequences
    
        align1 = ''
        align2 = ''
        i = len(a)
        j = len(b)
        
        while i > 0 or j > 0: # different from LCS to finish global alignment
            
            if Backtracks[i][j] == 0:            # pointer say go up a row
                align1 += a[i-1]        # add the base to align1  
                align2 += '-'           # add a gap to align2
                i -= 1                  # decrement row number
            
            elif Backtracks[i][j] == 1:          # pointer says up move, same as above for seq b
                align1 += '-'
                align2 += b[j-1]
                j -= 1
            
            elif Backtracks[i][j] == 2:          #diagonal move, respective bases to each align
                align1 += a[i-1]
                align2 += b[j-1]
                i -= 1; j -= 1
        
        #reverse path to have it start at first base for both seqs
        return align1[::-1], align2[::-1]
    #
    
    # Global Alignment Wrapper function   
    
    # Scoring parameters
    BLOSUM62 = Blosum62()
    sigma = 5
    
    Scores, backtrack = Scoring_and_Backtracks(a,b)  
    return Scores[-1][-1], list(Build_alignment(backtrack, a,b))

#####


with open("/Users/jasonmoggridge/Desktop/rosalind_test.txt",'r') as infile:
    a = infile.readline().strip()
    b = infile.readline().strip()

Score, align = global_alignment_BLOSUM62(a,b)

# Output is Score (overall for alignemnt)\n, SeqA aligned\n, seqB aligned

print(Score,'\n',align[0],'\n', align[1])
with open("/Users/jasonmoggridge/Desktop/rosalind_BA5D.txt",'w') as outfile:
    outfile.write(str(Score) + '\n')
    outfile.write(align[0] + '\n')
    outfile.write(align[1])
    outfile.close()
