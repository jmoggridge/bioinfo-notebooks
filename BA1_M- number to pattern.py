#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 03:31:04 2019

@author: jasonmoggridge

Implement NumberToPattern solved by 603
July 28, 2015, 9:03 p.m. by Rosalind Team
←→
Implement NumberToPattern
Convert an integer to its corresponding DNA string.

Given: Integers index and k.

Return: NumberToPattern(index, k).

Sample Dataset
45
4
Sample Output
AGTC
Extra Dataset

"""
def NumberToPattern(index, k):
    alpha = 'ACGT'
    string =''
    for deg in range(k-1,-1,-1): 
        string += alpha[index// 4**deg]    
        index = index % 4**deg
    return string
    




    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1m.txt', 'r')
index = int(f.readline().strip())
k = int(f.readline().strip())
    

index = 5437
k = 8
Pattern = NumberToPattern(index,k)    
