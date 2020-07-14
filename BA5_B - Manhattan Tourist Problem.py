#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 20:38:50 2020

@author: jasonmoggridge

Code Challenge: Find the length of a longest path in the
 Manhattan Tourist Problem.

Input: Integers n and m, followed by an n × (m + 1) matrix
 Down and an (n + 1) × m matrix Right. The two matrices are
 separated by the "-" symbol.

Output: The length of a longest path from source (0, 0) to
 sink (n, m) in the rectangular grid whose edges are defined 
 by the matrices Down and Right
"""

# Tourist fills in the matrix S of n+1,m+1 grid with longest path
# from 0,0 to that i,j. The longest path to source is stored in the
# last entry m of the last list n.

def tourist(n,m, Down, Right):
    
    # Source:    
    S[0][0] = 0
    
    # First Column:
    for i in range(1, n + 1):
        S[i][0] = S[i - 1][0] + Down[i - 1][0]

    # First row:
    for j in range(1, m + 1): 
        S[0][j] = S[0][j - 1] + Right[0][j - 1]
    
    # Subsequent rows, L to R filling S[i][j] of matrix:
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            North = S[i-1][j] + Down[i-1][j]
            West = S[i][j - 1] + Right[i][j-1]
            S[i][j] = max(North, West)
#   
   

#Parse input grids:
# n,m
# Down = n,m+1 array of weights for arrival from 'down' move
# - 
# right = n+1,m array of weight for arr from 'right' move


         
with open("/Users/jasonmoggridge/Desktop/rosalind_ba5b.txt",'r') as f:
    
    n,m = list(int(i) for i in f.readline().strip().split(' '))
    
    Down = list([] for i in range(n))
    for i in range(n):
        Down[i] = list(int(i) for i in f.readline().strip().split(' '))
    
    f.readline()
    
    Right = list([] for i in range(n+1))
    for i in range(n+1):
        Right[i] = list(int(i) for i in f.readline().strip().split(' '))    
    ##
    
    
# Initialize global var. S -matrix for longest paths to each node ij.
S = [[-1 for _ in range(m+1)]for _ in range(n+1)]    
# Call tourist to fill out grid with longest paths to i,j
tourist(n,m, Down, Right)
# Return longest path for sink n,m
longest_path = S[-1][-1]