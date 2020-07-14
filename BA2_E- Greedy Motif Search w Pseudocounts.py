#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 19:57:56 2019
@author: jasonmoggridge

BA2_E
    Implement Greedy Motif Search with Pseudocounts
    
    We encountered GreedyMotifSearch in “Implement GreedyMotifSearch”. 
    In this problem, we will power it up with pseudocounts.

Implement GreedyMotifSearch with Pseudocounts
    
    Given: 
        Integers k and t, followed by a collection of strings Dna.

    Return: 
        A collection of strings BestMotifs resulting from running 
        GreedyMotifSearch(Dna, k, t) with pseudocounts. If at any 
        step you find more than one Profile-most probable k-mer in
        a given string, use the one occurring first.
        """
###############################################################################


    
def Greedy_motif_search_w_pseudocounts(Seqs,k):
    """ Greedy motif search... 
        *select a kmer from first string
            *look at next string, find best k-mer using profile_matrix
            *update profile matrix*, move to next sequence
        *check motif alignment's total Hamming distance, save [Motifs set] 
        save [Motifs] with the smallest Hamming distance
            
    """
    
    def Build_profile_matrix(motifs, k):
        
        """Returns profile matrix (k,bases) given motifs"""
        ### PSEUDOCOUNTERS ENTERED WHILE PROFILE MATRIX IS CREATED    
        profile = [[1 for nt in range(len(alpha))] for position in range(k)]

        for p in range(k):
            for motif in motifs:
                profile[p][alpha.index(motif[p])] += 1
        return profile
        #
    
    
    def Profile_probable_motif(seq, profile):
    
        """Returns the most probable kmer to add to motifs, given profile matrix """
        
        kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
        best = -1
        
        for kmer in kmers:
            prob = 1
            for p in range(k):
                prob = prob * profile[p][alpha.index(kmer[p])]
            if prob > best:
                best = prob
                motif = kmer        
        return motif
        #
    
    
    def Mismatches(motifs):
    
        """ Returns the Hamming distance for  alignment of Motifs in a set from greedy"""
        
        consensus = [[0 for nt in range(4)] for position in range(k)]
        for p in range(k):
            for motif in motifs:
                consensus[p][alpha.index(motif[p])] += 1
                
        mismatches = 0
        for counts in consensus:
            for c in range(4): 
                if c != counts.index(max(counts)):
                    mismatches += counts[c]
        return mismatches
        #
        
        
    alpha ='ACGT'    
    best_mismatches = float('inf')
    Best_motifs = []
    
    for i in range(len(Seqs[0])-k+1):    
        Motifs = [Seqs[0][i:i+k]]
        for seq in Seqs[1:]:
            profile = Build_profile_matrix(Motifs, k)
            Motif = Profile_probable_motif(seq, profile)
            Motifs.append(Motif)    
        
        mismatches = Mismatches(Motifs)
        if  mismatches < best_mismatches:
            Best_motifs = Motifs
            best_mismatches = mismatches

    return Best_motifs


###############################################################################
  
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba2e.txt', 'r')
k,t = (int(i) for i in f.readline().strip().split(' ')); del(t)
Seqs = list(str(l.strip('\n')) for l in f.readlines())

Motifs = Greedy_motif_search_w_pseudocounts(Seqs, k)    
for motif in Motifs:
    print(motif)
    