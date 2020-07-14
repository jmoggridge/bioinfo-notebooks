#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 00:28:37 2020

@author: jasonmoggridge

Course 3: Edit distance problem

In 1966, Vladimir Levenshtein introduced the notion of the edit distance between two strings as the minimum number of edit operations needed to transform one string into another. Here, an edit operation is the insertion, deletion, or substitution of a single symbol. For example, TGCATAT can be transformed into ATCCGAT with five edit operations, implying that the edit distance between these strings is at most 5.

Edit Distance Problem
Find the edit distance between two strings.

Given: Two amino acid strings.

Return: The edit distance between these strings.

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



with open("/Users/jasonmoggridge/Desktop/rosalind_ba5g.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()


S, Pointers= Edit_Alignment(v,w)
edit_distance = -S[-1][-1]