#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 09:18:41 2020
@author: jasonmoggridge


CODE CHALLENGE: Solve the Peptide Identification Problem.
     Given: A space-delimited spectral vector Spectrum' and an amino acid string Proteome.
     Return: A substring of Proteome with maximum score against Spectrum'.

Sample Input:
    0 0 0 4 -2 -3 -1 -7 6 5 3 2 1 9 3 -8 0 3 1 2 1 8
    XZZXZXXXZXZZXZXXZ
Sample Output:
    ZXZXX

> peptide_identification.txt
KLEAARSCFSTRNE

dataset_11866_2
SOLUTION:
>   AATDVGALMYSPRW
"""

# -- Peptide Identification: Enumerate from proteome & Score by dot product
    # -- main function

def PeptideIdenfication(Proteome, Spectrum):

    mass = len(Spectrum) - 1
    low = mass // max(mass_aa.keys())
    hi =  mass // min(mass_aa.keys())
    best_score = -float('inf')
    
    for i in range(len(Proteome)-low+1):
        for j in range(i+low, i+hi):
            if j < len(Spectrum):
                if mass == Mass_Peptide(Proteome[i:j]):
                    score  = Dot_Product(Vector(Proteome[i:j]), Spectrum)
                    if score > best_score:
                        best_score = score
                        best_peptide = Proteome[i:j]
    
    return best_peptide


# -- Subroutines --
    
def Mass_Table():  
    mass_aa = {}; aa_mass = {}
    with open('data/integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_aa[int(mass)] = str(aa)
            aa_mass[str(aa)] = int(mass)
    return mass_aa, aa_mass


def Mass_Peptide(peptide):
    return sum(aa_mass[aa] for aa in peptide)


def Vector(peptide):
    prefixes = [sum(aa_mass[aa] for aa in peptide[:i]) for i in range(1,len(peptide)+1)]
    vector = [0] + [0 for _ in range(max(prefixes))]
    for prefix in prefixes:
        vector[prefix] = 1
    return vector


def Dot_Product(peptide, spectrum):    
    score = 0    
    for i in range(len(peptide)):
        score += peptide[i]*spectrum[i] 
    return score


# MAIN --
    
with open('data/rosalind_ba11f_2.txt','r') as infile:
    
    # Add a source node zero to the spectrum for (mass = 0)
    Spectrum = [0]+[int(i) for i in infile.readline().strip().split()]
    # Concatenated string of all peptides in proteome
    Proteome = infile.readline().strip()

mass_aa, aa_mass = Mass_Table()
peptide = PeptideIdenfication(Proteome, Spectrum)
print(peptide)




    