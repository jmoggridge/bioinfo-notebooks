#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 00:20:27 2019

@author: jasonmoggridge

BA3_A:
    
    Generate the k-mer Composition of a String 
    
    Given a string Text, its k-mer composition Compositionk(Text) 
    is the collection of all k-mer substrings of Text (including
    repeated k-mers). For example,

    Composition3(TATGGGGTGC) = {ATG, GGG, GGG, GGT, GTG, TAT, TGC, TGG}
    
    Note that we have listed k-mers in lexicographic order (i.e., how 
    they would appear in a dictionary) rather than in the order of their
    appearance in TATGGGGTGC. We have done this because the correct 
    ordering of the reads is unknown when they are generated.

String Composition Problem
    Generate the k-mer composition of a string.

    Given: An integer k and a string Text.

    Return: Compositionk(Text) (the k-mers can be provided in any order).
"""

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba3a.txt', 'r')
k = int(f.readline().strip())
s = str(f.readline().strip())

kmers=[]
for i in range(len(s)-k+1):
    kmers.append(s[i:i+k])
    
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_out.txt', 'w')
for kmer in kmers:
    f.write(str(kmer)+'\n')
f.close()
