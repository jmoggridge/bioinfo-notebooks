#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 13:32:27 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3e/

BA3_E:
    Construct the De Bruijn Graph of a Collection of k-mers

    Given an arbitrary collection of k-mers Patterns 
    (where some k-mers may appear multiple times), 
    we define CompositionGraph(Patterns) as a graph with 
    |Patterns| isolated edges. Every edge is labeled by a
    k-mer from Patterns, and the starting and ending nodes
    of an edge are labeled by the prefix and suffix of the
    k-mer labeling that edge. We then define the de Bruijn
    graph of Patterns, denoted DeBruijn(Patterns), by gluing
    identically labeled nodes in CompositionGraph(Patterns),
    which yields the following algorithm.

    DEBRUIJN(Patterns):
        represent every k-mer in Patterns as an isolated edge
        between its prefix and suffix
        glue all nodes with identical labels, yielding the 
        graph DeBruijn(Patterns)
        
        return DeBruijn(Patterns)


De Bruijn Graph from k-mers Problem:
    
    Construct the de Bruijn graph from a collection of k-mers.
    
    Given: A collection of k-mers Patterns.
    
    Return: The de Bruijn graph DeBruijn(Patterns), in the form of an adjacency list.    
    

"""
def DeBruijn(kmers):
    Adj = {}
    for kmer in kmers:
        Adj[kmer[:-1]] = []
    for kmer in kmers:
        Adj[kmer[:-1]].append(kmer[1:])
    return Adj
###
    

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba3e.txt', 'r')
kmers = list(str(l.strip('\n')) for l in f.readlines())

Adj =DeBruijn(kmers)
A = []
for adj in Adj:
    adjs = ','.join(Adj[adj])
    A.append(str(adj + ' -> ' + adjs))
    
wr = open('//Users/jasonmoggridge/Desktop/rosalind_outputBA3_E.txt', 'w')
for line in A:
    wr.write(line+'\n')
wr.close()