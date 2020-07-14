#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 19:24:08 2020
@author: jasonmoggridge



The following algorithm constructs ColoredEdges(P) for a genome P.
 In this pseudocode, we will assume that an n-element array (a1, . . . , an)
 has an invisible (n + 1)-th element that is equal to its first element, i.e., an+1 = a1.
 
ColoredEdges(P)
     Edges ← an empty set
     for each chromosome Chromosome in P
          Nodes ← ChromosomeToCycle(Chromosome)
          for j ← 1 to |Chromosome|
               add the edge (Nodes2j, Nodes2j +1) to Edges
     return Edges

Code Challenge: Implement ColoredEdges.
Input: A genome P.
Output: The collection of colored edges in the genome graph of P in the form (x, y).


Sample Input:

(+1 -2 -3)(+4 +5 -6)
Sample Output:

(2, 4), (3, 6), (5, 1), (8, 9), (10, 12), (11, 7)
"""



    

def colored_edgelist(Chromosomes):        
    
    def Chromosome_to_Cycle(Chromosome):
    
        Nodes =[False for _ in range(len(Chromosome)*2)]
        for j in range(1,len(Chromosome)+1):
            if Chromosome[j-1] > 0:
                Nodes[2*j-2] = 2*Chromosome[j-1] - 1
                Nodes[2*j-1] = 2*Chromosome[j-1]
            else:
                Nodes[2*j-2] = -2*Chromosome[j-1] 
                Nodes[2*j-1] = -2*Chromosome[j-1] - 1
        return Nodes
    #    
       
    col_edges = []
    for chromosome in Chromosomes:
        Nodes = Chromosome_to_Cycle(chromosome)
        for i in range(1,len(Nodes)-1,2):
            col_edges.append((Nodes[i], Nodes[i+1]))
        col_edges.append((Nodes[-1],Nodes[0]))
    print(', '.join(str(i) for i in col_edges))
    return col_edges



with open("/Users/jasonmoggridge/Desktop/dataset_8222_7.txt",'r') as infile:
    line = infile.readline().strip('(').strip(')')
    line = line[:-2]
    Chromosomes = [[int(j) for j in i.split(' ')] for i in line.split(')(')]
    infile.close(); del (line) 
    
edgelist = colored_edgelist(Chromosomes)