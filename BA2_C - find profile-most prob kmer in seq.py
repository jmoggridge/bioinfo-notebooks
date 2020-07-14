#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 16:13:30 2019

@author: jasonmoggridge

BA2_C

Find a Profile-most Probable k-mer in a String :
    
    Given a profile matrix P, 
    we can evaluate the probability of every k-mer in a string S
    and find a Profile-most probable k-mer in S, 
    i.e., a k-mer that was most likely to have been generated
    by P among all k-mers in S.
    
    For example, ACGGGGATTACC is the Profile-most probable
    12-mer in GGTACGGGGATTACCT. Indeed, every other 12-mer 
    in this string has probability 0.

    In general, if there are multiple Profile-most probable k-mers in Text,
    then we select the first such k-mer occurring in Text.

Profile-most Probable k-mer Problem

    Find a Profile-most probable k-mer in a string.

    Given: A string Text, an integer k, and a 4 × k matrix Profile.

    Return: A Profile-most probable k-mer in Text. (If multiple answers exist, you may return any one.)

Sample Dataset
    ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT
    5
    0.2 0.2 0.3 0.2 0.3
    0.4 0.3 0.1 0.5 0.1
    0.3 0.3 0.5 0.2 0.4
    0.1 0.2 0.1 0.1 0.2
Sample Output
    CCGAG

"""


def profile_most_pro_kmer(S,k,P):
    """looks at all k-mers in text and returns their p(seq|model P)
    p = 1* (each probability in path across matrix for k positions
    """
    alpha = 'ACGT'                          # need alpha to index base to matrix row: P[alpha.index(base)]
    kmer_prob ={}                           # dict to keep track of best results
    
    for i in range(len(S)-k+1):             # i indexes the kmer
        
        kmer = S[i:i+k]                     # 
        if kmer not in kmer_prob:           # add kmer to dict, set to p=zero
            kmer_prob[S[i:i+k]] = 0         
        score = 1                           # need to start score with p=1
        
        for j in range(0,k):
            #j iterates over kmer 
            base = alpha.index(kmer[j])     # index bases
            score = score * P[base][j]      # multiplies by P matrix(base, pos) for each base,pos
        
        if score > kmer_prob[kmer]:         # update if better score
            kmer_prob[kmer] = score
            print(kmer, score)    
    
    p_max_kmers = []                        # get best prob kmer in dict
    for kmer in kmer_prob.keys():
        if kmer_prob[kmer] == max(kmer_prob.values()):
            p_max_kmers.append(kmer)
            
    return p_max_kmers
##
    
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba2c.txt', 'r')
S = str(f.readline().strip())
k = int(f.readline().strip())
P = [[float(i) for i in l.strip('\n').split(' ')] for l in f.readlines()]

pmax = profile_most_pro_kmer(S,k,P)
for p in pmax:
    print('\n\n',p)
