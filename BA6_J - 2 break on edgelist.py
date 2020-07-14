#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 01:04:22 2020
@author: jasonmoggridge

1.6 CS: Solving the 2-Break Sorting Problem

The following pseudocode describes how 
2-Break(i1 , i2 , i3 , i4 ) transforms a genome graph.



2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4)
     remove colored edges (i1, i2) and (i3, i4) from GenomeGraph
     add colored edges (i1, i3) and (i2, i4) to GenomeGraph
     return GenomeGraph
     
Code Challenge: Implement 2-BreakOnGenomeGraph.

Input: 
    The colored edges of a genome graph GenomeGraph, 
    followed by indices i1 , i2 , i3 , and i4 .

Output: 
    The colored edges of the genome graph resulting 
    from applying the 2-break operation 
    2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 ).

Sample Input:

(2, 4), (3, 8), (7, 5), (6, 1)
1, 6, 3, 8

Sample Output:

(2, 4), (3, 1), (7, 5), (6, 8)
"""

with open("/Users/jasonmoggridge/Desktop/rosalind_ba6j.txt",'r') as infile:
    edges = infile.readline().strip()
    edges = edges[1:-1].split('), (')
    edges = [tuple(int(x) for x in i.split(', ')) for i in edges]
    
    indices = infile.readline().strip()
    v1,u1,v2,u2 = [int(i) for i in indices.split(', ')]
            
def two_break( edges, v1,u1, v2,u2 ):
    if (v1, u1) in edges:
        edges.remove((v1, u1))
    else:
        edges.remove((u1, v1))
    if (v2, u2) in edges:
        edges.remove((v2, u2))
    else:
        edges.remove((u2, v2))
    
    edges.append((v1, v2))
    edges.append((u1, u2))
    print(', '.join(str(i) for i in edges))
    return edges

edges = two_break(edges,v1,u1,v2,u2)    



