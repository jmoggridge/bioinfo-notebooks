#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 20:56:21 2019

@author: jasonmoggridge

Hammond Distance between two Sequences:
    We say that position i in k-mers p1 … pk and q1 … qk
    is a mismatch if pi ≠ qi. 
    For example, CGAAT and CGGAC have two mismatches. 
    The number of mismatches between strings p and q is
    called the Hamming distance between these strings a
    nd is denoted HammingDistance(p, q).

Hamming Distance Problem
Compute the Hamming distance between two DNA strings.

Given: Two DNA strings.

Return: An integer value representing the Hamming distance.
"""

def HammingDistance(p, q):
    H=0
    for i in range(len(p)):
        if p[i]!=q[i]:
            H +=1
    return H


f = open('//Users/jasonmoggridge/Desktop/dataset_9_3.txt', 'r')
p,q = str(f.readline().strip()),str(f.readline().strip())

H=0
for i in range(len(p)):
    if p[i]!=q[i]:
        H +=1
        
print(H)