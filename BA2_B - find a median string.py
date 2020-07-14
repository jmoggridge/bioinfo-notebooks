#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:18:17 2019

@author: jasonmoggridge

BA2_B: Median string


Find a Median String :
    
    Given a k-mer Pattern and a longer string Text,
    we use d(Pattern, Text) to denote the minimum 
    Hamming distance between Pattern and any k-mer in Text,

        d(Pattern,Text)=min(HammingDistance(Pattern,Pattern′)) for all slices

    Given a k-mer Pattern and a set of strings Dna = {Dna1, … , Dnat}, 
    we define d(Pattern, Dna) as the sum of distances between Pattern and 
    all strings in Dna,
        
     basically as above except taking the sum of distances over all strings
     
        d(Pattern,Dna)=∑(i,t) d(Pattern,Dna[i]).


    Our goal is to find a k-mer Pattern that minimizes d(Pattern, Dna) 
    over all k-mers Pattern, the same task that the Equivalent Motif 
    Finding Problem is trying to achieve. We call such a k-mer a median
    string for Dna.

Median String Problem
    Find a median string.

    Given: An integer k and a collection of strings Dna.

    Return: A k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern. (If multiple answers exist, you may return any one.)

Sample Dataset
k=    3
DNA={    
    AAATTGACGCAT
    GACGACCACGTT
    CGTCAGCGCCTG
    GCTGAGCACCGG
    AGTACGGGACAG
    }

Sample Output
    GAC
"""

def median_string(Seqs,k):
    
    # memo-ized Hammond distance between two eq len seqs
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
            ##
   
    # lists all DNA permutations of len k
    def all_kmers(k):
        import itertools
        alpha = "ACGT"
        all_kmers = [''.join(x) for x in itertools.product(alpha, repeat=k)]
        return all_kmers
        ##
    
    memo_H ={}
    kmers = list(all_kmers(k))
    median = ('',float('inf'))
    for kmer in kmers:
        score = 0
        for seq in Seqs:
            score += min((HammondD(kmer,seq[i:i+k])) for i in range(len(seq)-k+1))                     
        if score <= median[1]:
            median = (kmer, score)
    
    return median
    ##
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba2b.txt', 'r')
k = int(f.readline().strip())
Seqs = list(l.strip('\n') for l in f.readlines())
median = median_string(Seqs,k)
print(median)

