#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 17:09:17 2019
@author: jasonmoggridge

BA_4F:
    Compute the Score of a Cyclic Peptide Against a Spectrum

To generalize the Cyclopeptide Sequencing Problem from “Find a Cyclic Peptide with Theoretical Spectrum Matching an Ideal Spectrum” to handle noisy spectra, we need to relax the requirement that a candidate peptide’s theoretical spectrum must match the experimental spectrum exactly, and instead incorporate a scoring function that will select the peptide whose theoretical spectrum matches the given experimental spectrum the most closely. Given a cyclic peptide Peptide and a spectrum Spectrum, we define Score(Peptide, Spectrum) as the number of masses shared between Cyclospectrum(Peptide) and Spectrum. Recalling our example above, if

>Spectrum = {0, 99, 113, 114, 128, 227, 257, 299, 355, 356, 370, 371, 484},
then Score(NQEL, Spectrum) = 11.

The scoring function should take into account the multiplicities of shared masses, i.e., how many times they occur in each spectrum. For example, suppose that Spectrum is the theoretical spectrum of NQEL; for this spectrum, mass 242 has multiplicity 2. If 242 has multiplicity 1 in the theoretical spectrum of Peptide, then 242 contributes 1 to Score(Peptide, Spectrum). If 242 has larger multiplicity in the theoretical spectrum of Peptide, then 242 contributes 2 to Score(Peptide, Spectrum).

Cyclic Peptide Scoring Problem
Compute the score of a cyclic peptide against a spectrum.

Given: An amino acid string Peptide and a collection of integers Spectrum.

Return: The score of Peptide against Spectrum, Score(Peptide, Spectrum).

"""
    
mass = {
    'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,'I':113,'L':113,\
    'N':114,'D':115,'E':129,'K':128,'Q':128, 'M':131,'H':137,'F':147,\
    'R': 156,'Y': 163,'W': 186\
    }
##    


def cyclospectrum(protein):

    spec = [ 0, int( sum( mass[aa] for aa in protein ))]
    cycle = protein*2

    for i in range(len(protein)):
        for j in range(i+1, i+len(protein)):
            spec.append( sum( mass[aa] for aa in cycle[i:j] ))
            
    return spec
###

def score(peptide, spectrum):
    score = 0    
    for fragment in cyclospectrum(peptide):        
        if fragment in spectrum:
            score +=1
            spectrum.remove(fragment) 
    return score

#    
#f = open("/Users/jasonmoggridge/Desktop/rosalind_ba4f.txt", 'r')
#peptide = str(f.readline().strip())
#experimental = [int(i) for i in f.readline().split(' ')]

peptide = 'MAMA'
experimental = [0, 57, 71, 71, 71, 104, 131, 202, 202, 202, 256, 333, 333, 403, 404]
score = score(peptide, experimental)