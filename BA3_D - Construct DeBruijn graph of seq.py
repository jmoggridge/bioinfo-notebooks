#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 13:08:27 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3d/

BA3-D:
    Construct the De Bruijn Graph of a String
    
    
    Given a string Text, Path_Graph_k(Text) 
    is the path consisting of |Text| - k + 1 edges,
    where the i-th edge of this path is labeled by
    the i-th k-mer in Text and the i-th node of the
    path is labeled by the i-th (k - 1)-mer in Text.
    
    The de Bruijn graph DeBruijn_k(Text) is formed by
    gluing identically labeled nodes in PathGraphk(Text).

De Bruijn Graph from a String Problem:
    
    Construct the de Bruijn graph of a string.

    Given: 
        An integer k and a string Text.

    Return:
        DeBruijn_k(Text), in the form of an adjacency list.
    
"""
def DeBruijnize(seq,k):
    kmers = []
    for i in range(len(seq)-k+1):
        kmers.append(seq[i:i+k])
    
    Adj = {}
    for kmer in kmers:
        Adj[kmer[:-1]]=[]
    for kmer in kmers:
        Adj[kmer[:-1]].append(kmer[1:])
    
    return Adj


f = open('//Users/jasonmoggridge/Desktop/dataset_199_6.txt', 'r')
k = int(f.readline().strip())
seq = str(f.readline().strip())
Adj = DeBruijnize(seq,k)

A =[]
for adj in Adj:
    adjs = ','.join(Adj[adj])
    A.append(str(adj + ' -> ' + adjs))
    
    
wr = open('//Users/jasonmoggridge/Desktop/rosalind_outputBA3_D.txt', 'w')
for line in A:
    wr.write(line+'\n')
wr.close()