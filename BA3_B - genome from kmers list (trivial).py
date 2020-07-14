#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 12:41:20 2019

@author: jasonmoggridge

BA3_C : reconstruct genome path from list of kmers (trivial)
"""

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba3b.txt', 'r')
Kmers = list(str(l.strip('\n')) for l in f.readlines())

seq = Kmers[0]
for Kmer in Kmers[1:]:
    seq += Kmer[-1]
print(seq)