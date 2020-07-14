#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 23:21:28 2019

@author: jasonmoggridge


Find the Most Frequent Words with Mismatches in a String 

←→:
    We defined a mismatch in “Compute the Hamming Distance
    Between Two Strings”. We now generalize “Find the Most 
    Frequent Words in a String” to incorporate mismatches 
    as well.

    Given strings Text and Pattern as well as an integer d,
    we define Countd(Text, Pattern) as the total number of \
    occurrences of Pattern in Text with at most d mismatches.
    For example, Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4
    because AAAAA appears four times in this string with at 
    most one mismatch: AACAA, ATAAA, AAACA, and AAAGA. 
    Note that two of these occurrences overlap.
    
    A most frequent k-mer with up to d mismatches in Text is
     simply a string Pattern maximizing Countd(Text, Pattern)
     among all k-mers. Note that Pattern does not need to actually
     appear as a substring of Text; for example, AAAAA is the most
     frequent 5-mer with 1 mismatch in AACAAGCTGATAAACATTTAAAGAG,
     even though AAAAA does not appear exactly in this string. 
     Keep this in mind while solving the following problem.

Frequent Words with Mismatches Problem
        Find the most frequent k-mers with mismatches in a string.
        
        Given: A string Text as well as integers k and d.
        
        Return: All most frequent k-mers with up to d mismatches in Text.

Sample Dataset
1)
    ACGTTGCATGTCGCATGATGCATGAGAGCT
    4 1

Sample Output
    GATG ATGC ATGT 
    
2)
    Input:
CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC
10 2
    Output:
GCACACAGAC GCGCACACAC
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

def most_freq_kmers_w_mismatches(genome,k, threshold):

    # Get a set of all kmers from the input string.
    kmers = []
    for i in range(len(genome)-k+1):
        kmers.append(genome[i: i+k])

    # Generate all possible kmers, of len k
    import itertools
    all_kmers = [''.join(x) for x in itertools.product("ACGT", repeat=k)]
    
    #Initialize a dict and count how often kmer is within H <= d to a kmer from text
    counts = dict(zip(all_kmers,(0 for i in all_kmers)))
    for kmer in kmers:
        for a_kmer in all_kmers:
            if HammondD(kmer, a_kmer) <= threshold:
                counts[a_kmer]+=1
    
    # Get key of most frequently matched kmer
    most_freq=[]
    for kmer in all_kmers:
        if counts[kmer] == max(counts.values()):
            most_freq.append(kmer)            

    return most_freq
###
   
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1i.txt', 'r')
genome = str(f.readline().strip())
(k, threshold) = (int(i) for i in f.readline().strip().split(' '))
memo_H = {}
most_freq = most_freq_kmers_w_mismatches(genome, k, threshold)

string =''
for f in most_freq:
    string += str(f) + ' '
print(string)
#    
#


