#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 22:06:32 2019
@author: jasonmoggridge


BA2G - GIBBS SAMPLER:
    
    Works similar to random motif search except only replaces one motif in set
    during each loop.
    When selecting motif from randomly selected seq[i] in Seqs,
    new motif is sampled from seq[i]'s kmer based on p-distribution of
    motifs, given profile matrix from all seqs except seq[i]
    
    So making the profile matrix uses all seqs except the one being edited
    The new motif is chosen based on the probability given profile
    iterating over inner loop N times, instead of until reaches local maxima
        -has ability to get past poorly-scoring stretches in solution-space,
        unlike greedy.
    Scores are based on consensus alignment hamming distance. Low is better.

"""

import random


"""Generates a set of motifs by selecting a random kmer from each sequence"""        

def random_kmers(Seqs,k):
    motifs = [] 
    for seq in Seqs:
        x = random.randint(0, len(seq)-k)
        motifs.append(seq[x:x+k])
    return motifs
    ##


"""Returns profile matrix (k,bases) w pseudocounters given motifs"""

def Build_profile_matrix(motifs, k):
    
    # adds pseudocount of 1 to each (pos, base) entry in matrix    
    profile = [[1/(len(motifs)+4) for nt in range(len(alpha))] for position in range(k)]
    for p in range(k):
        for motif in motifs:
            profile[p][alpha.index(motif[p])] += 1/(len(motifs)+4)    
    return profile
    ##


"""Returns a random motif from seq based on profile probability weights"""
def profile_rand_motif(seq,k,profile):
    
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    c_weights = []
    c=0
    
    # generate cumulative probability array for weighted 'roulette' to select motif
    for kmer in kmers:
        prob = 1
        for p in range(k):
            prob = prob * profile[p][alpha.index(kmer[p])] 
        c+=prob
        c_weights.append(c)
    
    # j  is the random variable to land in roulette bins of array, if smaller than cum. value, take
    j = random.uniform(0,c_weights[-1])
    c_weights.append(float('inf'))
    for i in range(len(c_weights)):
        if c_weights[i] >= j:   
            return kmers[i]
   ##         
    

""" Returns the Hamming distance for  alignment of Motifs in a set from greedy"""
    
def Hamming(motifs):

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
    ##
    
  
    
###############################################################################
    
"""Wrapper function for Gibbs sampler motif searcher"""    
    
def Gibbs_Sampler(Seqs, k, N):        

    Best = float('inf')
    Best_motifs = []
    Motifs = random_kmers(Seqs,k)

    for _ in range(N):

        i = random.randint(0,len(Seqs)-1)
        profile_these_motifs = Motifs[:i] + Motifs[i+1:]
        Profile = Build_profile_matrix(profile_these_motifs, k)    
        Motifs[i] = profile_rand_motif(Seqs[i],k, Profile)    
        score = Hamming(Motifs)

        if score < Best:
            Best_motifs = Motifs
            Best = score

    return Best_motifs, Best
    
###############################################################################

    
f = open('//Users/jasonmoggridge/Desktop/dataset_163_4.txt', 'r')
k, x, N = (int(i) for i in f.readline().strip().split(' ')); del(x)
Seqs = list(str(l.strip('\n')) for l in f.readlines())
alpha ='ACGT'    

best_score = float('inf')
best_motifs =[]
for x in range(25):
    G_motifs, score = Gibbs_Sampler(Seqs, k, N)
    if score < best_score:
        print('best score', score, G_motifs)
        best_score = score
        best_motifs = G_motifs
        
for motif in best_motifs:
    print(motif)