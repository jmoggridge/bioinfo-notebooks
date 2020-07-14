#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 00:06:12 2019
@author: jasonmoggridge

BA3-I:
    Find a k-Universal Circular String:

        A k-universal circular string is a circular string 
        that contains every possible k-mer constructed over
        a given alphabet.

k-Universal Circular String Problem:
    
    Find a k-universal circular binary string.

    Given: An integer k.
    
    Return: A k-universal circular string. 
    (If multiple answers exist, you may return any one.)
###    
ALGO:
    
    Generate a set of all kmers over (alpha, k)
    Construct a DeBruijn graph for all (k-1)-mers from this set
    Use a dfs to find the Eulerian cycle guaranteed to exist
    Decode string from circuit by discarding last (k-1) edges in circuit
    
    
        
"""



import itertools
import sys; sys.setrecursionlimit(10**9) 

###

def DeBruijn_graph(kmers):        
    """Returns a DeBruijn graph for a set of overlapping kmers"""    
    graph = {}
    for kmer in kmers:
        if kmer[:-1] not in graph:
            graph[kmer[:-1]] = [kmer[1:]]
        else:
            graph[kmer[:-1]].append(kmer[1:])
    return graph
#
    

def dfs(v, graph, circuit, stack):
    """Eularian path/circuit dfs"""
    
    if len(graph[v]) == 0:
        circuit.append(v)
        if len(stack) > 0:
            dfs(stack.pop(), graph, circuit, stack)
    else:
        stack.append(v)
        dfs(graph[v].pop(-1), graph, circuit, stack)
    ##

def kmers_to_seq(Kmers):
    """Convert path of kmers back to a sequence"""    
    seq = Kmers[0]
    for Kmer in Kmers[1:]:
        seq += Kmer[-1]
    return seq
    ##
##
    
# k is the length of the strings the must be contained as subseqs of DeBruijn seq
k= 4
#alpha is binary
alpha = ('0','1')

#create all kmers from the alpha and length of words
all_kmers = [''.join(str(b) for b in kmer) for kmer in itertools.product(alpha,repeat=k)]
debruijn = DeBruijn_graph(all_kmers)

#start at the all zeroes k-1mer
start = str('0'*(k-1))

# dfs to give the eulerian circuit as a list
circuit = [] ; stack = []
dfs(start, debruijn, circuit, stack)
circuit.reverse()

#decoding the DeBruijn sequence from the eulerian circuit

s = circuit[0]
for c in circuit[:-(k-1)]:
    s += str(c[-1])
    

print(s+s[:k-2])

#
#
#binary_list = []
#for i in range(16):
#binary_list.append('{0:04b}'.format(i))