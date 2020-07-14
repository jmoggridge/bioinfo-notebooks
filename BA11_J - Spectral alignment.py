#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 03:49:26 2020

@author: jasonmoggridge

working dir:
    /Users/jasonmoggridge/Dropbox/Rosalind/Coursera_textbook_track/Course4_Phylogeny
    
    
    solves sample
    solves extra sample
    fails test
    
"""
#
## toy dataset masses
#aa_mass = {'X':4,'Z':5}
#mass_aa = {4:'X', 5:'Z'}
#
#


def getMass_Tables():  
    
    #look-up tables for (AA:mass)
    
    mass_aa = {}; aa_mass = {}
    with open('data/integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_aa[int(mass)] = str(aa)
            aa_mass[str(aa)] = int(mass)
    return mass_aa, aa_mass

mass_aa, aa_mass = getMass_Tables()
#

def Spectral_Alignment(peptide, spectrum, k):
    
    
        
        
        # compute the diff array for amino acids in peptide
    # compute mass array for cumulative mass -> row headers in matrix
    diff = {0:0}
    masses = [0]
    for aa in peptide:
        masses.append(masses[-1] + aa_mass[aa])
        diff[masses[-1]] = aa_mass[aa]
        
    # init the score and prev arrays for dynamic programming/memo'ing longest path in DAG
    score = {}
    prev = {}
    for i in masses:
        for j in range(len(spectrum)):
            for t in range(0,k+1):
                score[(i,j,t)] = -float('inf')
    score[(0,0,0)] = 0
    
    
    
    
    # t is the number of edits in current path, allow up to k edits..
    for t in range(0, k+1):
        
    #    print('\n\n\n for number of edits t:', t)
    
        #iterate over i row -> AAi in peptide (indexed by mass of substring pept[0:i])
        for i in masses[1:]:
        
    #        print('\n\t For sub peptide i = ', i)
            
            for j in range(1, len(spectrum)):
        
                # predecessor in same t-level of graph; ie. default to the diagonal edge
                if j-diff[i] >= 0:
                    score[(i,j,t)] = spectrum[j] + score[(i-diff[i], j - diff[i], t)]
                    prev[(i,j,t)] = (i-diff[i], j - diff[i], t)
    
                else:
                    score[(i,j,t)] = -float('inf')
                    prev[(i,j,t)] = None
        #        
                # if not on bottom t-level, j predecessors on t-1 level to consider
                if t > 0:
                    for jj in range(0,j):
                        if jj != j-diff[i]:
    
                            s = spectrum[j] + score[(i-diff[i], jj, t-1)]
                            if s > score[(i,j,t)]:
                                score[(i,j,t)] = s
                                prev[(i,j,t)] = (i-diff[i], jj, t-1)
                        
    #            if score[(i,j,t)] > -float('inf'):
    #                print('\ti,j,t',(i,j,t),'score[(i,j,t)] ->', score[(i,j,t)],'\tprev[(i,j,t)] ->', prev[(i,j,t)])
        
    
    # Find the best score in the bottom/right of matrix, take best of k levels
    best = - float('inf')
    for edits in range(0,k+1):
        if score[(masses[-1], len(spectrum)-1, edits)] > best:
            (i,j,t) = (masses[-1], len(spectrum)-1, edits)
            best = score[(masses[-1], len(spectrum)-1, edits)] 
    
    # build sequence in reverse using backpointers
    seq = ''
    while prev[(i,j,t)]:
        
        (pi,pj,pt) = prev[(i,j,t)]
        aa = mass_aa[i - pi]
        
        # if AA mass was modified -> identify mass difference, add to AA string
        actual = j - pj   
        if (i-pi) != actual:
            delta = actual - (i-pi)
            if delta >0:
                aa += '(+' + str(delta) + ')'
            else:
                aa += '(' + str(delta) + ')'
    
        # tacks new aa on to the front of the seq strubg
        seq = aa + seq
    
        # break loop once have reached the source
        (i,j,t) = (pi,pj,pt)
        if (i,j,t) == (0,0,0):
            break
        
    print('peptide ->\n\t',seq)


    
    
    
    
    
#main

with open("data/dataset_11866_14.txt",'r') as file:
    peptide = file.readline().strip()
    spectrum = [0]+ [int(i) for i in file.readline().strip().split()]
    k = int(file.readline().strip())
    
seq = Spectral_Alignment(peptide, spectrum, k)
























