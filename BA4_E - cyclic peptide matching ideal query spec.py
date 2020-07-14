#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 17:35:55 2019
@author: jasonmoggridge


In “Compute the Number of Peptides of Given Total Mass”, we first encountered the problem of reconstructing a cyclic peptide from its theoretical spectrum; this problem is called the Cyclopeptide Sequencing Problem and is given below. It is solved by the following algorithm.

    CYCLOPEPTIDESEQUENCING(Spectrum)
        Peptides ← a set containing only the empty peptide
        while Peptides is nonempty
            Peptides ← Expand(Peptides)
            for each peptide Peptide in Peptides
                if Mass(Peptide) = ParentMass(Spectrum)
                    if Cyclospectrum(Peptide) = Spectrum
                        output Peptide
                    remove Peptide from Peptides
                else if Peptide is not consistent with Spectrum
                    remove Peptide from Peptides
                    
                    
Cyclopeptide Sequencing Problem:
    
    Given an ideal experimental spectrum, find a cyclic peptide
    whose theoretical spectrum matches the experimental spectrum.
    
    Given: 
        A collection of (possibly repeated) integers Spectrum
        corresponding to an ideal experimental spectrum.
    
    Return: 
        Every amino acid string Peptide such that 
        Cyclospectrum(Peptide) = Spectrum
        (if such a string exists).
    
sample:
    spec = '0 113 128 186 241 299 314 427'
Sample Output
   cyclopeptides -> 186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
"""
#import sys
#sys.setrecursionlimit(10**8) 

mass = {
        'G':57,'A':71,'S':87,'P':97,'V':99,'T':101,'C':103,\
        'I':113,'N':114,'D':115,'E':129,'K':128, 'M':131,'H':137,\
        'F':147,'R': 156,'Y': 163,'W': 186
        }
##

def cyclopeptide_sequencing(experimental):
    
    """ Finds all cyclopeptides that would generate experimental spectral (idealized)
    
    Algo:
      Expand peptides:  
          Progressively generates larger possible linear subpeptides
          starting with [0] mass, try all possible aa's

      Linear & Cyclo- Spectrum fxs: 
          return a list of fragment weights
          for ideal mass spec from query peptide
          
      Bound:
          Bound expanded subpeptides list by checking whether
          peptide's linear spectrum is consistent with experimental
          
          Once mass of peptides has reached largest frag in experimental:
            need to check the cyclospectrum of the peptide
                if it matches experimental, save this solution
        
        """        
    
    def expand_peptides(peptides):

        new_peptides = []
        for peptide in peptides:
            for aa in mass:
                new_peptides.append(peptide + aa)

        return new_peptides
        ##
    
    def linear_spectrum(peptide):
    
        spectrum =[0]
        for i in range(len(peptide)):
            for j in range(i+1, len(peptide)+1):
                spectrum.append(sum(mass[aa] for aa in peptide[i:j]))
        return spectrum
        ##
        
    def cyclo_spectrum(peptide):
    
        spectrum = [0, int(sum(mass[aa] for aa in peptide))]
    
        cycle = peptide*2
        for i in range(len(peptide)):
            for j in range(i+1, i+len(peptide)):
                spectrum.append(sum(mass[aa] for aa in cycle[i:j]))
            
        return sorted(spectrum)
        ##    
        
    def bound(peptides):
            
        bound_peptides = []
        for peptide in peptides:
            
            consistent = True
            peaks = list(experimental)
            
            spec = linear_spectrum(peptide)
            for fragment in spec:
                if fragment not in peaks:
                    consistent = False
                    peaks = list(experimental)
                    break
                else:
                    peaks.remove(fragment)        
            if consistent:
                bound_peptides.append(peptide)
                
                if max(spec) == sorted(experimental)[-1]:
                    if cyclo_spectrum(peptide) == experimental:
                        print('found a matching cyclopeptide', peptide)

                        final_peptides.append(peptide)
    
        return bound_peptides
        ##
    
    final_peptides=[]
    peptides = ['']
    while peptides:
        peptides = expand_peptides(peptides)
        peptides = bound(peptides)
    
    return final_peptides

##

"""DATA INPUT AND OUTPUT"""


f = open("/Users/jasonmoggridge/Desktop/dataset_100_6.txt",'r')
experimental = [int(i) for i in f.readline().split(' ')]

cyclopeptides = cyclopeptide_sequencing(experimental)
 
###
    
f = open("/Users/jasonmoggridge/Desktop/rosalind_ba4e_output.txt",'w')

for peptide in cyclopeptides:
    masses = [str(mass[aa]) for aa in peptide]
    string = '-'.join(masses)
    f.write(string+' ')
f.close()