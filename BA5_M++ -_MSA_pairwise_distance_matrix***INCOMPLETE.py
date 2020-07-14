
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 17:52:34 2020
@author: jasonmoggridge

Return the best scoring pair-wise alignment from 
a set of input sequences.

"""

Seqs = []
with open("/Users/jasonmoggridge/Desktop/rosalind_ba5m.txt",'r') as infile:
    for seq in infile.readlines():
        Seqs.append('-'+ str(seq.strip()))
    infile.close()
    
sigma = 1
    

import numpy as np

def Pairwise_Alignment_Score( v, w, sigma, scoring):
        
    S = np.zeroes((len(w), len(v)), dtype=np.int16)
    P = -1 * np.ones((len(w), len(v)), dtype=np.int16)
    
#    for i in range(1, len(v)):
#        S[i][0] = -i*sigma         # need to penalize all the gaps in first col
#        Backtracks[i][0] = 0
#        
#    for i in range(1, len(w)):       # penalize all gaps in first row
#        Scores[0][i] = -i * sigma
#        Backtracks[0][i] = 1
#        
#        for i in range(1, len(a) +1):
#            for j in range(1, len(b) +1):
#    
#                match = BLOSUM62[(a[i-1], b[j-1])]
#                
#                score = max(                # pick longest path from 3 options
#                 Scores[i -1][j] - sigma,
#                 Scores[i][j-1] - sigma,
#                 Scores[i-1][j-1] + match)
#                Scores[i][j] = score
#                
#                if score == Scores[i][j-1] - sigma:
#                    Backtracks[i][j] = 1
#                elif score == Scores[i-1][j] - sigma:
#                    Backtracks[i][j] = 0
#                elif score == Scores[i-1][j-1] + match:
#                    Backtracks[i][j] = 2
#                    
#        return Scores, Backtracks
#     