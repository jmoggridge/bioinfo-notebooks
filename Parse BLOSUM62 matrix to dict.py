#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 19:00:30 2020

@author: jasonmoggridge


Parse BLOSUM62
"""

#Blosum62
with open("/Users/jasonmoggridge/Desktop/BLOSUM62.txt", 'r') as file:
    AA = file.readline().strip().split('  ')
    arrays = []
    for _ in AA:
        array = list(file.readline().strip().split(' '))
        while '' in array:
            array.remove('')
        array = list(int(i) for i in array[1:])
        arrays.append(array)
    del(array)

BLOSUM62 = {}
for i in range(len(AA)):
    for j in range(len(AA)):
        BLOSUM62[(AA[i],AA[j])] = arrays[i][j]
        
        
        
        
        
def PAM250():
    
    with open("/Users/jasonmoggridge/Desktop/PAM250.txt", 'r') as file:
        AA = file.readline().strip().split('  ')
        arrays = []
        for _ in AA:
            array = list(file.readline().strip().split(' '))
            while '' in array:
                array.remove('')
            array = list(int(i) for i in array[1:])
            arrays.append(array)
    
    pam250 = {}
    for i in range(len(AA)):
        for j in range(len(AA)):
            pam250[(AA[i],AA[j])] = arrays[i][j]
            
    return pam250