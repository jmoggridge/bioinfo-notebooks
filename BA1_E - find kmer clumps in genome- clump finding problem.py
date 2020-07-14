#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:57:20 2019

@author: jasonmoggridge

Find Patterns Forming Clumps in a String 

←→:
    Given integers L and t, a string Pattern forms an (L, t)-clump
    inside a (larger) string Genome if there is an interval of 
    Genome of length L in which Pattern appears at least t times.
    
    For example, TGCA forms a (25,3)-clump in the following 
    Genome: gatcagcataagggtcccTGCAATGCATGACAAGCCTGCAgttgttttac.

Clump Finding Problem:
Find patterns forming clumps in a string.

Given: A string Genome, and integers k, L, and t.

Return: All distinct k-mers forming (L, t)-clumps in Genome.

Sample Dataset

CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC
5 75 4

Sample Output
CGACA GAAGA AATGT

Extra Datasets


Takes more than a couple seconds for a long string, could maybe speed it up by searching for clumps more 
efficiently through data structure and string properties, skipping indices where kmer cant appear..

"""

    
def find_clumps_in_subseq(genome, k, L, t):
    print('Working on finding:\n\t(L,t)-clumps of: k=',k,'-mers\n\t-Appearing at least t=',t,'times\n\t-in subseq of len L=',L)
    #generates a list of distinct kmers in genome
    
    def generate_distinct_kmers(genome,k):
    
        kmers = []
        for i in range(len(genome)-k+1):
            kmer = genome[i: i+k]
            if kmer not in kmers:
                kmers.append(kmer)
        return kmers
    ##

    kmers = generate_distinct_kmers(genome, k)
    clumped = []

    for i in range(len(genome) - L + 1):                # for all subseqs of len L in g
        
        subseq = genome[i:i+L]                          # this is the subseq slice
        counts = dict(zip(kmers,(0 for i in kmers)))    # init counts for all kmers

        for i in range(L-k+1):                          # slide across subseq, count the kmer there
            counts[subseq[i:i+k]] +=1                   # up the count of kmer found
        
        for kmer in counts:
            if counts[kmer] >= t:                       # check to see which k-mers
                if kmer not in clumped:                 # have count > t threshold
                    clumped.append(kmer)                # if unique, add to list

    s = ''
    for kmer in clumped:
        s += str(kmer) + ' '
    
    print('\n\n And here are your clumps sire:\n\n', s)
    return s
##
            
            
f = open('//Users/jasonmoggridge/Desktop/dataset_4_5.txt', 'r')
genome = str(f.readline().strip())
(k,L,t) = (int(i) for i in f.readline().strip().split(' '))
clumped = find_clumps_in_subseq(genome, k, L, t)    
    
### FASTER ALGORITHM

def ClumpFinding(genome,k,L,t):
        
    # Create a dictionary for {kmer:freq}
    
    Kmer_count = {}
    for i in range(len(genome)-k+1):
        kmer = genome[i: i+k]
        if kmer not in Kmer_count:
            Kmer_count[kmer] = 0
    
    #Generate initial freq array and set first done and new kmers to update
    
    Clumps = []
    
    for i in range(L-k+1):
        kmer = genome[i:i+k]
        Kmer_count[kmer] += 1
        if Kmer_count[kmer] >= t:
            if kmer not in Clumps:       
                Clumps.append(kmer)
                print(kmer, 'is a clump in first window')
    
    #iterate over all steps in L sliding window, update just outgoing/incoming kmer
    for i in range(1,len(genome) - L+1):
        
        outgoing = genome[i-1: i+k-1]         
        Kmer_count[outgoing]-= 1
        
        incoming = genome[i+L-k : i+L]        
        Kmer_count[incoming] += 1
        if Kmer_count[incoming] >= t:
            if incoming not in Clumps:       
                Clumps.append(incoming) 
                print(incoming, 'is a new found clump')
     
    return Clumps

            
    
with open("/Users/jasonmoggridge/Dropbox/Rosalind/textbook_track/Course1_Hidden_messages_in_DNA/e_coli.txt",'r') as f:
    genome = str(f.readline().strip())

k = 9       # k = 9-mers, occurring...
L = 500     # within a L = 500 bp window
t = 3       # at least t = 3 times
 
clumps = ClumpFinding(genome,k,L,t)