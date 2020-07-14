#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 00:16:50 2019

@author: jasonmoggridge

LINEAR SPECTRUM OF A PEPTIDE FUNCTION


"""

def linear_spectrum(peptide):

    spectrum =[0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide)+1):
            spectrum.append(sum(mass[aa] for aa in peptide[i:j]))
    return spectrum
    ##
    
mass = {
    'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,\
    'N':114,'D':115,'E':129,'K':128,'Q':128, 'M':131,'H':137,'F':147,\
    'R': 156,'Y': 163,'W': 186\
    }
##

peptide ='SEQIWEGAQNIIPRKLYWCTNSFMQKEHPEVPAHFP'

spec = linear_spectrum(peptide)
string =''
for mass in sorted(spec):
    string += str(mass)+' '
    
print(string)