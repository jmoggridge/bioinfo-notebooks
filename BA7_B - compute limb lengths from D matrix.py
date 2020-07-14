#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 22:32:43 2020
@author: jasonmoggridge

Limb Length Problem
Find the limb length for a leaf in a tree.


We now have an algorithm for solving the Limb Length Problem. 
For each j, we can compute LimbLength(j) by finding the 
minimum value of (D_{i,j} + D_{j,k} - D_{i,k}) / 2
over all pairs of leaves i and k.

Given: An integer n, followed by an integer j between 0 and n - 1,
 followed by a space-separated additive distance matrix D (whose elements are integers).

Return: The limb length of the leaf in Tree(D) corresponding to row
 j of this distance matrix (use 0-based indexing).

Sample Dataset
4
1
0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0
"""
# Limb len for node j in distance matrix D using triangle inequality
def Limb(D,j):
    
    n = len(D)
    limb_len = float('inf')
    for i in range(n):
        if j != i:
            for k in range(n)[i+1:]:
                if k != j:
                    limb_candidate = (D[i][j] + D[j][k] - D[i][k]) /2
                    if limb_candidate < limb_len:
                        limb_len = int(limb_candidate)

    return limb_len
    
# Limb len for entire set of leafs in distance matrix
    
def Limb_Lengths(D):
    
    n = len(D)
    limb_len = [float('inf') for _ in range(n)]
    
    for j in range(n):
        for i in range(n-1):
            if j != i:
                for k in range(n)[i+1:]:
                    if k != j:
                        limb_candidate = (D[i][j] + D[j][k] - D[i][k]) /2
                        if limb_candidate < limb_len[j]:
                            limb_len[j] = int(limb_candidate)

    return limb_len

# Main

with open("/Users/jasonmoggridge/Desktop/rosalind_test2.txt",'r') as infile:
    n = int(infile.readline())
    j = int(infile.readline())

    D = []
    lines = infile.readlines()
    for line in lines:
        D.append(list(int(i) for i in line.strip().split()))
    del(lines, line)
    infile.close()
    