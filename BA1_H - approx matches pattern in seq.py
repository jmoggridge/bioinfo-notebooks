#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 23:09:40 2019

@author: jasonmoggridge


Find All Approximate Occurrences of a Pattern in a String 

←→:
    We say that a k-mer Pattern appears as a substring of Text
    with at most d mismatches if there is some k-mer substring
    Pattern' of Text having d or fewer mismatches with Pattern,
    i.e., HammingDistance(Pattern, Pattern') ≤ d. 
    
    Our observation that a DnaA box may appear with slight variations
    leads to the following generalization of the Pattern Matching Problem.

Approximate Pattern Matching Problem:
    
    Find all approximate occurrences of a pattern in a string.
    
    Given: Strings Pattern and Text along with an integer d.
    
    Return: All starting positions where Pattern appears as a substring of Text with at most d mismatches.

Sample Dataset
ATTCTGGA
CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC
3
Sample Output
6 7 26 27 78


Basic pattern matching whenever d <= threshold

"""

def HammondD(p,q):
    H=0
    for i in range(len(p)):
        if p[i]!=q[i]:
            H +=1 
    return H
    ##
    
def Approx_matches_in_seq(genome, pattern, threshold):

    start_pos = []

    for i in range(len(genome)-len(pattern)+1):
        if HammondD(pattern, genome[i:i+len(pattern)]) <= threshold:
            start_pos.append(i)


    return start_pos

    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1h.txt', 'r')
pattern = str(f.readline().strip())
genome = str(f.readline().strip())
d = int(f.readline().strip())
approx_matches =Approx_matches_in_seq(genome, pattern, d)
count = len(approx_matches)
