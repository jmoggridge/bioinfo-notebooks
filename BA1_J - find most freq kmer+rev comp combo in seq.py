#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 20:14:30 2019
@author: jasonmoggridge

BA1_J:
    
Find Frequent Words with Mismatches and Reverse Complements:
    
    We now extend “Find the Most Frequent Words with Mismatches in a String”
    to find frequent words with both mismatches and reverse complements.
    Recall that RC refers to the reverse complement of Pattern.


Frequent Words with Mismatches and Reverse Complements Problem

    Find the most frequent k-mers (with mismatches and reverse complements) in a DNA string.
    
    Given: A DNA string Text as well as integers k and d.
    
    Return: All k-mers Pattern maximizing the sum Countd(Text, Pattern) 
    + Countd(Text, RC) over all possible k-mers.

Sample Dataset
    ACGTTGCATGTCGCATGATGCATGAGAGCT
    4 1
Sample Output
    ATGT ACAT

"""
    

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
    
def RC(kmer):
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    kmer_rc = ''
    for bp in kmer:
        kmer_rc = complement[bp] + kmer_rc
    return kmer_rc
    ##



def most_freq_kmers_w_mismatches(genome,k, threshold):

    # Get a set of all kmers from the input string.
    kmers = []
    for i in range(len(genome)-k+1):
        kmers+= [genome[i: i+k], RC(genome[i: i+k])]

    # Generate all possible kmers, of len k
    import itertools
    all_kmers = [''.join(x) for x in itertools.product("ACGT", repeat=k)]
    
    #Initialize a dict and count how often kmer is within H <= d to a kmer from text
    counts = dict(zip(all_kmers,(0 for i in all_kmers)))

    for kmer in kmers:
        for a_kmer in all_kmers:
            if HammondD(kmer, a_kmer) <= threshold:
                counts[a_kmer]+=1
    
    # Get key of most frequently matched kmer/RC combo in string
    
    most_freq = []
    target = max(counts[kmer] + counts[RC(kmer)] for kmer in all_kmers)
    
    for kmer in counts.keys():
        if counts[kmer] + counts[RC(kmer)] == target:
            most_freq.append(kmer)            

    return most_freq
    ###


memo_H = {}    

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1j.txt', 'r')
genome = str(f.readline().strip())
(k, threshold) = (int(i) for i in f.readline().strip().split(' '))

most_freq = most_freq_kmers_w_mismatches(genome, k, threshold)

string =''
for f in most_freq:
    string += str(f) + ' '
print(string)
#    
#

