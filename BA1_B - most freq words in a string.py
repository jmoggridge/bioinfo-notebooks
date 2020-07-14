#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 18:55:52 2019

@author: jasonmoggridge



Find the Most Frequent Words in a String solved by 1308
July 28, 2015, 9:11 p.m. by Rosalind Team
←→
We say that Pattern is a most frequent k-mer in Text if it maximizes Count(Text, Pattern) among all k-mers. For example, "ACTAT" is a most frequent 5-mer in "ACAACTATGCATCACTATCGGGAACTATCCT", and "ATA" is a most frequent 3-mer of "CGATATATCCATAG".

Frequent Words Problem
Find the most frequent k-mers in a string.

Given: A DNA string Text and an integer k.

Return: All most frequent k-mers in Text (in any order).

Sample Dataset
"""

#text = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
#k = 4
#most_freq_kmers = ['CATG', 'GCAT']


def most_freq_kmers(text, k):
    
    kmers = {}
    
    for i in range(len(text)-k+1):
        kmer = text[i: i+k]
        if kmer not in kmers:
            kmers[kmer] = 1
        else:
            kmers[kmer] += 1
                    
    freq = [] 
    for kmer in kmers.keys():
        if kmers[kmer] == max(kmers.values()):
            freq.append(kmer)

    return freq
#   #        


f = open("/Users/jasonmoggridge/Desktop/rosalind_ba1b.txt",'r')
text = f.readline().strip()
k = int(f.readline().strip())
f.close()

freq = most_freq_kmers(text, k)

output = open("/Users/jasonmoggridge/Desktop/rosalind_ba1b_out.txt",'w')
for kmer in freq:
    print(str(kmer))
    output.write(str(kmer)+' ')
    
output.close()        
