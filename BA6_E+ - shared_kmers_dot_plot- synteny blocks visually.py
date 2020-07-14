#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 19:49:01 2020

@author: jasonmoggridge

Find All Shared k-mers of a Pair of Strings 

We say that a k-mer is shared by two genomes if either the
 k-mer or its reverse complement appears in each genome. 
 In Figure 1 are four pairs of 3-mers that are shared by
 "AAACTCATC" and "TTTCAAATC".

A shared k-mer can be represented by an ordered pair (x, y),
 where x is the starting position of the k-mer in the first
 genome and y is the starting position of the k-mer in the
 second genome. For the genomes "AAACTCATC" and "TTTCAAATC",
 these shared k-mers are (0,4), (0,0), (4,2), and (6,6).

Shared k-mers Problem
Given two strings, find all their shared k-mers.

Given: An integer k and two strings.

Return: All k-mers shared by these strings, in the form 
of ordered pairs (x, y) corresponding to starting positions
 of these k-mers in the respective strings.

Sample Dataset

3
AAACTCATC
TTTCAAATC


Sample Output

(0, 4)
(0, 0)
(4, 2)
(6, 6)
"""


def Shared_Kmers(a, b): 
    # a,b are DNA sequences.Returns a list of indices x,y for location in a,b

    def Reverse_Complement(seq):
        rc = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
        return str(''.join(rc[base] for base in seq))[::-1]
        
    kmers = {}
    for i in range(len(a)-k+1):    
        kmer = a[i:i+k] 
        if kmer not in kmers.keys():
            kmers[kmer] = [i]
            kmers[Reverse_Complement(kmer)] = [i]
        else:
            kmers[kmer].append(i)
            kmers[Reverse_Complement(kmer)].append(i)
            
    shared_kmers = [] 
    for y in range(len(b)-k+1):
        kmer = b[y:y+k]
        if kmer in kmers.keys():
            for x in kmers[kmer]:
                shared_kmers.append((x,y))

    return list(set(shared_kmers))


# MAIN
    

with open("/Users/jasonmoggridge/Desktop/dataset_289_5.txt",'r') as infile:
    k = int(infile.readline().strip())
    a = infile.readline().strip()
    b = infile.readline().strip()
    infile.close()

with open("/Users/jasonmoggridge/Desktop/outfile.txt",'w') as outfile:
    for kmer in Shared_Kmers(a,b):
        outfile.write(str(kmer)+'\n')
    outfile.close
    
    
    
