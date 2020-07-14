#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:01:35 2020
@author: jasonmoggridge

BA6_F - Chromosome to Cycle

The following pseudocode bypasses the intermediate step of assigning
 “head” and “tail” nodes in order to transform a single circular
 chromosome 
 Chromosome = (Chromosome1, . . . , Chromosome_n) into 
 a cycle represented as a sequence of integers 
 
 Nodes = (Nodes1, . . . , Nodes2n).

ChromosomeToCycle(Chromosome)
     for j ← 1 to |Chromosome|
          i ← Chromosomej
          if i > 0
               Node2j−1 ←2i−1
               Node2j ← 2i
          else
               Node2j−1 ← -2i
               Node2j ←-2i−1
     return Nodes
Chromosome To Cycle Problem
Solve the Chromosome To Cycle Problem.

Given: A chromosome Chromosome containing n synteny blocks.

Return: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.

Sample Dataset
(+1 -2 -3 +4)
Sample Output
(1 2 4 3 6 5 7 8)

"""


with open("/Users/jasonmoggridge/Desktop/rosalind_ba6f.txt",'r') as infile:
    line = infile.readline().strip('(').strip(')')
    Chromosome = [int(i) for i in line.split(' ')]
    infile.close()

def Chromosome_to_Cycle(Chromosome):
    
    Nodes =[False for _ in range(len(Chromosome)*2)]
    for j in range(1,len(Chromosome)+1):
        if Chromosome[j-1] > 0:
            Nodes[2*j-2] = 2*Chromosome[j-1] - 1
            Nodes[2*j-1] = 2*Chromosome[j-1]
        else:
            Nodes[2*j-2] = -2*Chromosome[j-1] 
            Nodes[2*j-1] = -2*Chromosome[j-1] - 1
    print(' '.join(str(i) for i in Nodes))
    return Nodes

Nodes = Chromosome_to_Cycle(Chromosome)