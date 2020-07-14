#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 08:26:45 2020
@author: jasonmoggridge

BA11_E - > peptide sequencing problem


    Recall 'the change problem' -> dynamic programming for optimal path in DAG.
    
            #def DPChange(money, coins):
            #    
            #    minCoins = {0:0}
            #    for m in range(1,money+1):
            #        minCoins[m] = float('inf')
            #        for coin in coins:
            #            if m >= coin:
            #                if minCoins[m-coin]+1 < minCoins[m]:
            #                    minCoins[m] = minCoins[m-coin]+1
            #    print(minCoins)
            #    return minCoins[money]
    ###
DP strategy -> peptide seq'n
    - Use Amino acids masses' as edges between fragments i,j, (aa instead of coins)
    - total spectrum mass instead of money (the final subproblem)
    
    - Memoize subproblem solutions (scores, backedge pointer) for each mass (same way)
    - iterating from mass =1 -> mass = mass_peptide (m= len(S)) (same)
    - need to add backtracking to rebuild peptide
    
sample data
     with open('data/peptide_sequencing.txt','r') as infile:
        Spectrum = [0]+[int(i) for i in infile.readline().strip().split()]
    #
    #Output
    #GGPGGPGGAGG

 toy masses
    #mass_aa = {4:'X',5:'Z'}
    #aa_mass ={'X':4,'Z':5}
    
   with open('data/dataset_11813_10.txt','r') as infile:
           > GPAAGLAGL
           
           
"""


# Peptide Sequencing function --

def DP_peptide_sequencing(Spectrum):
    
    # Builds mass->aa and aa->mass hashes
    def Mass_Table():      
    
        mass_aa = {}; aa_mass = {}    
        with open('data/integer_mass_table.txt', 'r') as file:
            for line in  file.readlines():
                [aa,mass] = line.strip().split()
                mass_aa[int(mass)] = str(aa)
                aa_mass[str(aa)] = int(mass)    
        return mass_aa, aa_mass
    ##
    
    # init mass tables, scores, pointers for DP
    mass_aa, aa_mass = Mass_Table() 
    m = len(Spectrum)
    scores = [0]+ [-float('inf') for _ in range(m-1)]
    pointers = [False for _ in range(m)]
    
    # DP -> memo-ize best paths to each mass in spectrum, save score and pointer
    for mass in range(1, m):
        for res in mass_aa.keys():
            if mass >= res:
                if scores[mass] < scores[mass - res] + Spectrum[mass]:
                    scores[mass] = scores[mass - res] + Spectrum[mass]
                    pointers[mass] = mass - res
    
    # DP -> backtrack pointers from sink -> source, gets backedge amino acid from ∆(mass j,i)
    # reconstruct peptide from end to start (make sure to reverse it!)
    j = m-1
    peptide = ''
    while j != 0:
        i = pointers[j]
        peptide += mass_aa[j-i]
        j = i 
        
    print(peptide[::-1])
#
    
# Main --
    
with open('data/rosalind_ba11e.txt','r') as infile:
    Spectrum = [0]+[int(i) for i in infile.readline().strip().split()]

DP_peptide_sequencing(Spectrum)

# Solution to rosalind dataset -> rosalind_ba11e.txt -> HGGAGGPCSAP





