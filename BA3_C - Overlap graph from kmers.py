#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 12:47:12 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3c/

BA3_C:

    Construct the Overlap Graph of a Collection of k-mers
    
    In this chapter, we use the terms prefix and suffix to refer
    to the first k − 1 nucleotides and last k − 1 nucleotides of
    a k-mer, respectively.

    Given an arbitrary collection of k-mers Patterns, we form a 
    graph having a node for each k-mer in Patterns and connect 
    k-mers Pattern and Pattern' by a directed edge if Suffix(Pattern)
    is equal to Prefix(Pattern'). The resulting graph is called 
    the overlap graph on these k-mers, denoted Overlap(Patterns).
    
    Overlap Graph Problem:
        
        Construct the overlap graph of a collection of k-mers.
            
        Given: A collection Patterns of k-mers.
        
        Return: The overlap graph Overlap(Patterns), in the form of an adjacency list.
"""


def Overlap_Graph(Kmers):    
    Adj = []
    for kmer1 in Kmers:
        for kmer2 in Kmers:
            if kmer1[1:] ==  kmer2[:-1]:
                Adj.append(str(kmer1 + ' -> ' + kmer2))
    return Adj
###

f = open('//Users/jasonmoggridge/Desktop/dataset_198_10.txt', 'r')
Adj = Overlap_Graph(list(str(l.strip('\n')) for l in f.readlines()))
wr = open('//Users/jasonmoggridge/Desktop/rosalind_outputBA3_C.txt', 'w')
for line in Adj:
    wr.write(line+'\n')
wr.close()