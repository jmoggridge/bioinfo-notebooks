#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 00:03:57 2019

@author: jasonmoggridge

BA2H - Implement Distance(Between Pattern And Strings)

    The first potential issue with implementing MedianString
    from “Find a Median String” is writing a function to 
    compute d(Pattern, Dna) = ∑ti=1 d(Pattern, Dnai), 
    the sum of distances between Pattern and each string 
    in Dna = {Dna1, ..., Dnat}. 
    This task is achieved by the following pseudocode.


    DistanceBetweenPatternAndStrings(Pattern, Dna)
        k ← length(Pattern)
        distance ← 0
        for each string Text in Dna:
            HammingDistance ← ∞
            for each k-mer Pattern’ in Text:
                if HammingDistance > HammingDistance(Pattern, Pattern’)
                    HammingDistance ← HammingDistance(Pattern, Pattern’)
            distance ← distance + HammingDistance
        return distance

Compute DistanceBetweenPatternAndStrings
Find the distance between a pattern and a set of strings.

Given: A DNA string Pattern and a collection of DNA strings Dna.

Return: DistanceBetweenPatternAndStrings(Pattern, Dna).
"""
######


def Distance(Pattern, Seqs):
    
    ##
    def Hamming(seq1, seq2): #distance for pattern -> kmer of equal length from seq
        distance = 0
        for i in range(len(seq1)):
            if seq1[i] != seq2[i]:
                distance += 1
        return distance
    ##

    k = len(Pattern)
    distance_pattern = 0
    
    for dna in Seqs:
        hamming_d = float('inf')
        kmers = [dna[i:i+k] for i in range(len(dna)-k+1)]
        
        for kmer in kmers:
            score = Hamming(kmer, Pattern)
            if score < hamming_d:
                hamming_d = score
        distance_pattern += hamming_d
    
    return distance_pattern

###

    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba2h.txt', 'r')
Pattern = (str(f.readline().strip()))
DNAs = list(str(i) for i in f.readline().strip().split(' '))
###
D = Distance(Pattern, DNAs)