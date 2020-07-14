#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 01:40:04 2019
@author: jasonmoggridge

BA4_C:
    Spectrum of cyclic peptide.
    
    Generate all possible fragments from the cyclic peptide
    Sum amino acid masses for each fragment
    compile fragment weights into sorted list
    don't ignore duplicate weights 
    
    Generating Theoretical Spectrum Problem: Generate the theoretical spectrum of a cyclic peptide.
        
        Input: An amino acid string Peptide.
        Output: Cyclospectrum(Peptide).
    
"""

def cyclospectrum(protein):
    
    # monoisotopic mass table
    mass = {
        'A':71.03711, 'C':103.00919, 'D':115.02694,\
        'E':129.04259, 'F':147.06841, 'G':57.02146,\
        'H':137.05891, 'I':113.08406,'K': 128.09496,\
        'L': 113.08406,'M': 131.04049,'N': 114.04293,\
        'P': 97.05276,'Q': 128.05858,'R': 156.10111,\
        'S': 87.03203,'T': 101.04768,'V': 99.06841,\
        'W': 186.07931,'Y': 163.06333 }
    
    for aa in mass:
        mass[aa] = int(round(mass[aa]))
    ###
    
    fragments=[protein]
    cycle =protein*2
    for i in range(len(protein)):
        for j in range(i+1, i+len(protein)):
            fragments.append(cycle[i:j])
        
    spectrum =[]
    for fragment in fragments:
        spectrum.append(sum(mass[aa] for aa in fragment))
        
    string = '0'
    for fragment in sorted(spectrum):
        string += ' ' + str(fragment)
        
    print(string)
    return string, spectrum
###
    
f = open("/Users/jasonmoggridge/Desktop/rosalind_ba4c.txt",'r')
protein = str(f.readline().strip())
protein ='IDCCYQHHHSEEG'
string,spectrum = cyclospectrum(protein)

f.close()

