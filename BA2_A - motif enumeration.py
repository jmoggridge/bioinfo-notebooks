#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 04:51:16 2019

@author: jasonmoggridge

BA2_A:

Implement 'Motif Enumeration': -> find all(k,d)-motifs that match all strings
    
    Given a collection of strings Dna and an integer d, 
    a k-mer is a (k,d)-motif if it appears in every string
    from Dna with at most d mismatches. The following algorithm
    finds (k,d)-motifs.

    MOTIFENUMERATION(Dna, k, d)
        Patterns ← an empty set []
        for each k-mer Pattern in Dna
            for each k-mer Pattern’ differing from Pattern by at most d
              mismatches
                if Pattern' appears in each string from Dna with at most d
                mismatches
                    add Pattern' to Patterns
        remove duplicates from Patterns
        return Patterns


Implanted Motif Problem

    Implement MotifEnumeration (shown above) to find all (k, d)-motifs in a collection of strings.
    
    Given: Integers k and d, followed by a collection of strings Dna.
    
    Return: All (k, d)-motifs in Dna.

Sample Dataset
    3 1
    ATTTGGC
    TGCCTTA
    CGGTATC
    GAAAATT
Sample Output
    ATA ATT GTT TTT

"""
 

def motif_enumeration(dna,k,d):
    
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
    #
    # K -> list of lists of kmers in each seq.  
    K = []
    for seq in dna:
        kmers = []
        for i in range(len(seq)-k+1):
            kmers.append(seq[i:i+k])
        K.append(kmers)
    
    # Generating d-'neighborhood': all the kmers with hammondD no greater than d to any seq kmer
    neighbourhood = []
    for seq in K:
        for kmer in seq:
            for neigh in all_kmers(k):
                if HammondD(kmer, neigh) <= d:
                    if neigh not in neighbourhood:
                        neighbourhood.append(neigh)
   
    # take a list of uniques from union of kmers and their d-neighborhood
    kmers =[]
    for seq in K:
        kmers += seq
    queries = list(set(kmers + neighbourhood))
    
    # Patterns -> find the kmers that have occurences in all strings (with mismatches<= d)
    patterns = []
    for query in queries:
        found = [False for seq in dna]
        for i in range(len(K)):
            for kmer in K[i]:
                if HammondD(kmer, query) <= d:
                    found[i] = True
        if all(found):
            patterns.append(query)
    return patterns

    ###
f = open('/Users/jasonmoggridge/Dropbox/dataset_156_8.txt', 'r')
k,d = (int(i) for i in f.readline().strip().split(' '))
DNA = [l.strip('\n') for l in f.readlines()]
memo_H ={}
patterns = motif_enumeration(DNA,k,d)            


f = open('//Users/jasonmoggridge/Desktop/rosalind_op.txt', 'w')
for p in patterns:
    f.write(p +' ')
f.close()
