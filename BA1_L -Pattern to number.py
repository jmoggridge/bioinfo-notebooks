#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 03:04:41 2019

@author: jasonmoggridge

Implement PatternToNumber solved by 636
July 28, 2015, 9:07 p.m. by Rosalind Team
←→
Implement PatternToNumber
Convert a DNA string to a number.

Given: A DNA string Pattern.

Return: PatternToNumber(Pattern).

Sample Dataset
AGT
Sample Output
11

"""

    
f = open('//Users/jasonmoggridge/Desktop/rosalind_ba1l.txt', 'r')
pattern = str(f.readline().strip())

##

def PatternToNumber(Pattern):
    alpha ='ACGT'
    poly =len(Pattern)-1
    number=0
    
    for i in range(0, len(Pattern)):
        number += 4**poly *alpha.index(Pattern[i])    
        poly -=1
    return number

##
Pattern ='ATGCAA'
number = PatternToNumber(Pattern)

