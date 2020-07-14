#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 06:13:10 2020
@author: jasonmoggridge

CODE CHALLENGE: Solve the Peptide Sequencing Problem.
     Given: A space-delimited spectral vector Spectrum'.
     Return: An amino acid string with maximum score against Spectrum'. For masses
     with more than one amino acid, any choice may be used.

Note: When a spectral vector Spectrum' = s1 ... sm is given, it does not have
 a zero-th element; in your implementations, you should assume that s0 is equal to zero.

Sample Input:
0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8
Sample Output:
XZZXX
"""

#mass_aa = {4:'X', 5:'Z'}
#aa_mass = {'X':4, 'Z':5}



def Mass_Table():

    mass_aa = {}
    aa_mass = {}

    with open('data/integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_aa[int(mass)] = str(aa)
            aa_mass[str(aa)] = int(mass)

    return mass_aa, aa_mass

def Vector(peptide):
    prefixes = [sum(aa_mass[aa] for aa in peptide[:i]) for i in range(1,len(peptide)+1)]
    vector = [0 for _ in range(max(prefixes))]
    for prefix in prefixes:
        vector[prefix-1] = 1
    return vector

with open('data/rosalind_ba11c.txt', 'r') as infile:
    peptide = infile.readline().strip()
#    peptide = 'SMQKTNMTIYAYYGIIEHDVCKVYDQYKWCNQFNDFI'

mass_aa, aa_mass = Mass_Table()
vector = Vector(peptide)

print((' ').join(str(i) for i in vector))


