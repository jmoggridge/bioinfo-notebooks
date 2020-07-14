#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 22:31:57 2019
@author: jasonmoggridge

BA2_F

    RANDOM MOTIF SEARCH -> find best set of conserved motifs in a set of seqs

    outer loop ->repeat 1000x to converge on best solution:
        -generate a set of random motifs
        -best score = infinity
        loop:
            generate a motif Profile matrix, given current Motifs
            select the profile-most probable motif from each string
            check if score of this set of motifs is < best score
            if so - store as best motifs, repeat loop
            else - repeat outerloop unless done
        output best set of motifs (lowest hamming distance for consensus alignment)
"""
import random


"""Returns profile matrix (k,bases) w pseudocounters given motifs"""

def Build_profile_matrix(motifs, k):
        
    profile = [[1/(len(motifs)+4) for nt in range(len(alpha))] for position in range(k)]
    for p in range(k):
        for motif in motifs:
            profile[p][alpha.index(motif[p])] += 1/(len(motifs)+4)    
    return profile
    ##


"""Returns the most probable kmers from each seq, given profile matrix """
        
def Profile_probable_motifs(Seqs, profile):
    motifs = []
    for seq in Seqs:    
        kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
        best = -1    
        for kmer in kmers:
            prob = 1
            for p in range(k):
                prob = prob * profile[p][alpha.index(kmer[p])]
            if prob > best:
                best = prob
                motif = kmer        
        motifs.append(motif)
    return motifs
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
    #
    
def random_kmers(Seqs,k):
    motifs = [] 
    for seq in Seqs:
        x = random.randint(0, len(seq)-k)
        motifs.append(seq[x:x+k])
    return motifs
    ##
    
def Random_motifs_searcher(Seqs, k):        

    Best = float('inf')
    Best_motifs = []
    Motifs = random_kmers(Seqs,k)

    improving = True
    while improving: 
#        print('improving....')
        Profile = Build_profile_matrix(Motifs, k)
        Motifs = Profile_probable_motifs(Seqs, Profile)
        score = Hamming(Motifs)
        if score < Best:
#            print('new best...improving')
            Best_motifs = Motifs
            Best = score
        else:
#            print('returning best...')
            return Best_motifs, Best
    
def random_motifs_wrapper(Seqs, k, repeat):
    
    GOAT = float('inf')
    GOAT_motifs = []   

    for _ in range(repeat):
#        print('Doing another repeat')
        R_Motifs, score = Random_motifs_searcher(Seqs,k)
#        print(R_Motifs, score)
        if score < GOAT:
#            print('new score: =', score)
            GOAT_motifs = R_Motifs
            GOAT = score

#    print('Returning GOAT')

    return GOAT_motifs
###############################################################################

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba2f.txt', 'r')
k,x = (int(i) for i in f.readline().strip().split(' ')); del(x)
Seqs = list(str(l.strip('\n')) for l in f.readlines())
alpha ='ACGT'    

repeat = 1000
Best_rand_motifs = random_motifs_wrapper(Seqs, k, repeat)

for motif in Best_rand_motifs:
    print(motif)
