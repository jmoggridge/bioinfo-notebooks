#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 02:51:16 2019
@author: jasonmoggridge

Generate the Frequency Array of a String 
←→:
    Given an integer k, we define the frequency array of a string Text
    as an array of length 4^k, where the i-th element of the array holds
    the number of times that the i-th k-mer (in the lexicographic order)
    appears in Text 

Computing a Frequency Array
    Generate the frequency array of a DNA string.
    
    Given: A DNA string Text and an integer k.
    
    Return: The frequency array of k-mers in Text.

Sample Dataset
    ACGCGGCTCTGAAA
    2
Sample Output
    2 1 0 0 0 0 2 2 1 2 1 0 0 1 1 0


"""
import itertools
def freq_array(genome, k):
    
    # Get a set of all kmers from the input string.

    alpha = "ACGT"
    kmers = [''.join(x) for x in itertools.product(alpha, repeat=k)]
    
    counts = [0 for j in kmers]
    for i in range(len(genome)-k+1):
        counts[kmers.index(genome[i: i+k])] += 1
    return counts
    ##
    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1k.txt', 'r')
genome = str(f.readline().strip())
k = int(f.readline().strip())

counts = freq_array(genome, k)
s =''
for i in counts:
    s += str(i)+ ' '
print(s)

print('\nDONE')