#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 02:00:39 2019

@author: jasonmoggridge

Find a Position in a Genome Minimizing the Skew solved by 925
July 28, 2015, 9:03 p.m. by Rosalind Team
←→
Define the skew of a DNA string Genome, denoted Skew(Genome), as the difference between the total number of occurrences of 'G' and 'C' in Genome. Let Prefixi (Genome) denote the prefix (i.e., initial substring) of Genome of length i. For example, the values of Skew(Prefixi ("CATGGGCATCGGCCATACGCC")) are:

0 -1 -1 -1 0 1 2 1 1 1 0 1 2 1 0 0 0 0 -1 0 -1 -2

Minimum Skew Problem
Find a position in a genome minimizing the skew.

Given: A DNA string Genome.

Return: All integer(s) i minimizing Skew(Prefixi (Text)) over all values of i (from 0 to |Genome|).
"""


S =open('//Users/jasonmoggridge/Desktop/rosalind_ba1f.txt', 'r').read().strip()
    
def min_skew_pos(genome):
    G = 0
    C = 0
    skew=[]
    for s in S:
        if s == 'G':
            G +=1
        elif s == 'C':
            C +=1
        skew.append(G-C)
    
    ans = []
    target = min(skew)
    for s in range(len(skew)):
        if skew[s] == target:
           ans.append(s+1) # add one bc zero-index of list 
    return ans, target    

min_skew, skew = min_skew_pos(S)
