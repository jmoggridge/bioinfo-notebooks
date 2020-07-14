#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:45:40 2019
@author: jasonmoggridge

Find All Occurrences of a Pattern in a String 
←→:
    In this problem, we ask a simple question: how many times
    can one string occur as a substring of another? Recall from
    “Find the Most Frequent Words in a String” that different
    occurrences of a substring can overlap with each other. 
    For example, ATA occurs three times in CGATATATCCATAG.

Pattern Matching Problem:
    Find all occurrences of a pattern in a string.
    
    Given: Strings Pattern and Genome.
    Return: All starting positions in Genome where Pattern 
    appears as a substring. Use 0-based indexing.


Sample Dataset

ATAT

GATATATGCATATACTT

Sample Output

1 3 9

"""

f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1d.txt', 'r')
Pattern = str(f.readline().strip())
Genome = str(f.readline().strip())


def pattern_in_string(Pattern, Genome):
    index = []
    for i in range(len(Genome)-len(Pattern)+1):
        if Genome[i:i+len(Pattern)] == Pattern:
            index.append(i)
    return index

ind = pattern_in_string(Pattern, Genome)
ans =''
for i in ind:
    ans += str(i)+' '
print('\n\n\nIndices of occurences of pattern in string (start positions, zero-indexed:\n')
print(ans)