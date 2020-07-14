#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 17 02:34:39 2020
@author: jasonmoggridge

The colored edges in the breakpoint graph of P and Q are given
by ColoredEdges(P) together with ColoredEdges(Q). Note that some
 edges in these two sets may connect the same two nodes, which
 results in trivial cycles.

Although we are now ready to solve the 2-Break Distance Problem,
 we will later find it helpful to implement a function converting
 a genome graph back into a genome.


GraphToGenome(GenomeGraph)
     P ← an empty set of chromosomes
     for each cycle Nodes in GenomeGraph
          Nodes﻿ ← sequence of nodes in this cycle (starting from node 1)
          Chromosome ← CycleToChromosome(Nodes)
          add Chromosome to P
     return P


 Code Challenge: Implement GraphToGenome.

Input: The colored edges ColoredEdges of a genome graph.
Output: The genome P corresponding to this genome graph
"""



### BA6I - Edgelist graph back to (chromosomes) (ie genome =[]):

def Graph_To_Genome(graph):

    cycles = []
    cycle = []
    while graph:
        edge = graph.pop(0)
        if not cycle:
            if edge[0] % 2 == 0:
                end = edge[0] - 1
            else:
                end = edge[0] + 1
        cycle += [e for e in edge]
        if cycle[-1] == end:
            cycles.append([cycle[-1]]+cycle[:-1])
            cycle = []

    Genome = []
    for cycle in cycles:
        Genome.append(Cycle_to_Chromosome(cycle))
    return Genome
