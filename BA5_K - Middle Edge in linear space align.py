#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:17:52 2020

@author: jasonmoggridge

Find a middle edge in the alignment graph in linear space.

Given: Two amino acid strings.

Return: A middle edge in the alignment graph of these strings,
 where the optimal path is defined by the BLOSUM62 scoring matrix 
 and a linear indel penalty equal to 5. Return the middle edge in the 
 form “(i, j) (k, l)”, where (i, j) connects to (k, l).
 

A space-efficient seq algorithm

    - use two column structure to store (_,j)(_,j-1) values for longest path
        - split into two subproblems for tracing best alignment
    - identify middle node
    - divide and conquer

Code Challenge: Solve the Middle Edge in Linear Space Problem (for protein strings).

Input: Two amino acid strings.

Output: A middle edge in the alignment graph in the form "(i, j) (k, l)",
 where (i, j) connects to (k, l). To compute scores, use the BLOSUM62
 scoring matrix and a (linear) indel penalty equal to 5.
 
 Middle Edge in Linear Space Problem

 
 """



def Blosum62():
    B62 = {}
    with open("/Users/jasonmoggridge/Dropbox/Rosalind/textbook_track/Course3_Alignments/BLOSUM62.txt", 'r') as file:
        AA = file.readline().strip().split('  ')
        arrays = []
        for _ in AA:
            array = list(file.readline().strip().split(' '))
            while '' in array:
                array.remove('')
            array = list(int(i) for i in array[1:])
            arrays.append(array)
    for i in range(len(AA)):
        for j in range(len(AA)):
            B62[(AA[i],AA[j])] = arrays[i][j]
    return B62





def Mid_Edge(v,w):

    def Paths_to_mid(v,w):
    
        S = [[0,-sigma]] + list([-sigma*i,-float('inf')] for i in range(1,len(v)+1))    
        for j in range(1,len(w)+1):
            for i in range(len(v)+1):
    
                match = mu[(v[i-1], w[j-1])]
    
                S[i][j%2] = max(\
                 S[i-1][j%2] - sigma,\
                 S[i][(j-1)%2] - sigma,\
                 S[i-1][(j-1)%2] + match)
    
        return list(S[i][j%2] for i in range(len(v)+1))
    
    mid_j = len(w)//2

    From_Source = Paths_to_mid(v, w[:mid_j])
    To_Sink = Paths_to_mid(v[::-1], w[len(w):mid_j-1:-1])
    To_Sink.reverse()
    
    best = -float('inf')
    for i in range(len(v)+1):
        len_ipath = From_Source[i] + To_Sink[i]
        if len_ipath >= best:
            best = len_ipath
            mid_i = i

    score = max(-sigma, mu[(v[mid_i], w[mid_j])])

    if score == -sigma:
        return  (mid_i, mid_j),(mid_i, mid_j+1)
    else:
        return (mid_i, mid_j),(mid_i+1, mid_j+1)





with open("/Users/jasonmoggridge/Desktop/rosalind_ba5k.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()

mu = dict(Blosum62())
sigma = 5

v = 'PLEASANTLY'
w ='MEASNLY'
Edge =  Mid_Edge(v,w)
print(' '.join(str(e) for e in Edge))
