#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 13:35:44 2019

@author: jasonmoggridge
"""

"""DeBruijn Directed Graphs:"""

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
    