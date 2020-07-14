#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 18:45:06 2020
@author: jasonmoggridge

BA6_G - CYCLE TO CHROMOSOME


The process described in “Implement ChromosomeToCycle” is in fact invertible, as described by the following pseudocode.

CycleToChromosome(Nodes)
     for j ← 1 to |Nodes|/2
          if Node2j−1 < Node2j
               Chromosomej ← Node2j /2
          else
               Chromosomej ← −Node2j−1/2
     return Chromosome
Cycle To Chromosome Problem
Solve the Cycle to Chromosome Problem.

Given: A sequence Nodes of integers between 1 and 2n.

Return: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.

Sample Dataset
(1 2 4 3 6 5 7 8)
Sample Output
(+1 -2 -3 +4)
"""




def Cycle_to_Chromosome(Nodes):    
    Chromosome = [False for _ in range(0, len(Nodes),2)]
    for i in range( 0, len(Nodes),2):
        if Nodes[i] < Nodes[i+1]:
            Chromosome[i//2] = (Nodes[i]+1)//2
        else:
            Chromosome[i//2] = -int((Nodes[i]+1)//2)
    return Chromosome





with open("/Users/jasonmoggridge/Desktop/dataset_8222_5.txt",'r') as infile:
    line = infile.readline().strip('(').strip(')')
    Nodes = [int(i) for i in line.split(' ')]
    infile.close(); del (line)

string = Cycle_to_Chromosome(Nodes)
