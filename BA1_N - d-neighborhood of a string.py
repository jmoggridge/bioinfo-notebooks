#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 04:35:37 2019

@author: jasonmoggridge


Generate the d-Neighborhood of a String :
    
    The d-neighborhood Neighbors(Pattern, d)
    is the set of all k-mers whose Hamming
    distance from Pattern does not exceed d.

Generate the d-Neighborhood of a String:
    
    Find all the neighbors of a pattern.
    
    Given: A DNA string Pattern and an integer d.
    
    Return: The collection of strings Neighbors(Pattern, d).

Sample Dataset
    ACG
    1
Sample Output
    CCG
    TCG
    GCG
    AAG
    ATG
    AGG
    ACA
    ACC
    ACT
    ACG
"""

###
def HammondD(p,q):
    
    if (p,q) in memo_H.keys():
        return memo_H[(p,q)]
    elif (q,p) in memo_H.keys():
        return memo_H[(q,p)]
    else:
        H=0
        for i in range(len(p)):
            if p[i]!=q[i]:
                H +=1 
        memo_H[(p,q)] = H
        memo_H[(q,p)] = H
        return H
    
def generate_all_kmers(k):

    import itertools
    alpha = "ACGT"
    kmers = [''.join(x) for x in itertools.product(alpha, repeat=k)]

    return kmers
###
#    
#f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1n.txt', 'r')
#pattern = str(f.readline().strip())
#d = int(f.readline().strip())

#pattern ='CCAGTCAATG'
#d=1
    

memo_H ={}

d_neighbors=[]
for kmer in generate_all_kmers(len(pattern)):
    if HammondD(kmer, pattern) <= d:
        d_neighbors.append(kmer)


f = open('//Users/jasonmoggridge/Desktop/rosalind_op.txt', 'w')
for d in d_neighbors:
    f.write(d +'\n')
f.close()
