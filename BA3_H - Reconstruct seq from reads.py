#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 23:22:44 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3h/

BA3_H:
    Reconstruct a String from its k-mer Composition
    
    String Reconstruction Problem:
        
        Reconstruct a string from its k-mer composition.
        
        Given: An integer k followed by a list of k-mers Patterns.
        
        Return: A string Text with k-mer composition equal to Patterns.
         (If multiple answers exist, you may return any one.)

ALGO:
StringReconstruction(Patterns)
    dB ← DeBruijn(Patterns)
    path ← EulerianPath(dB)
    Text﻿ ← PathToGenome(path)
    return Text
        
"""


def string_reconstruction(kmers):
        
    import sys
    sys.setrecursionlimit(10**6)     
    
    def DeBruijn_graph(kmers):        
        """Returns a DeBruijn graph for a set of overlapping kmers"""    
        Adj = {}
        for kmer in kmers:
            if kmer[:-1] not in Adj:
                Adj[kmer[:-1]] = [kmer[1:]]
            else:
                Adj[kmer[:-1]].append(kmer[1:])
        return Adj
    ##
    

    def get_source(Adj):   
        """Find the source node to initiate dfs at:"""
        deg = {}
        for v in Adj:       
            if v not in deg.keys():
                deg[v] = -len(Adj[v])
            else:
                deg[v] -= len(Adj[v])   
            for u in Adj[v]:
                if u not in deg.keys():
                    deg[u] = 1
                else:
                    deg[u] += 1   
        for v in deg:
            if deg[v] == -1:
                source = v    
        return source
        ##
        
    def dfs(v, graph, path, stack):
        """Eularian path/circuit dfs"""
        if v not in graph.keys():
            path.append(v)
            if len(stack) == 0:
                return
            else:
                dfs(stack.pop(), graph, path, stack) 
        elif len(graph[v]) == 0:
            path.append(v)
            if len(stack) == 0:
                return
            else:
                dfs(stack.pop(), graph, path, stack)
        else:
            stack.append(v)
            dfs(graph[v].pop(0), graph, path, stack)
        ##
                
    def kmers_to_seq(Kmers):
        """Convert path of kmers back to a sequence"""    
        seq = Kmers[0]
        for Kmer in Kmers[1:]:
            seq += Kmer[-1]
        return seq
        ##
    
    # String Reconstruction wrapper function:    
    debruijn = DeBruijn_graph(kmers)
    source = get_source(debruijn)    
    path = []
    stack =[]
    dfs(source,debruijn, path, stack)
    path.reverse()
    seq = kmers_to_seq(path)
    return seq
    #
#


f = open('//Users/jasonmoggridge/Desktop/rosalind_ba3h.txt', 'r')
k = int(f.readline().strip())
kmers = list(str(l.strip('\n')) for l in f.readlines())
seq = string_reconstruction(kmers)
print(seq)
    

        

