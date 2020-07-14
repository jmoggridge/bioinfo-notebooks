#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 23:32:58 2020

@author: jasonmoggridge


PSM-search algorithm -> find best matching peptides given spectra, proteome

PSMSearch(SpectralVectors, Proteome, threshold).
    PSMSet ← an empty set
    for each vector Spectrum' in SpectralVectors
          Peptide ← PeptideIdentification(Spectrum', Proteome)
          if Score(Peptide, Spectrum) ≥ threshold
              add the PSM (Peptide, Spectrum') to PSMSet
    return PSMSet
    
#"""



def getMass_Tables():  
    
    #look-up tables for (AA:mass)
    
    mass_aa = {}; aa_mass = {}
    with open('data/integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_aa[int(mass)] = str(aa)
            aa_mass[str(aa)] = int(mass)
    return mass_aa, aa_mass

mass_aa, aa_mass = getMass_Tables()
#

def Score(peptide, spectrum):

    # Computes dot-product of peptide vector and spectrum

    prefixes = [sum(aa_mass[aa] for aa in peptide[:i]) for i in range(1,len(peptide)+1)]
    vector = [0 for _ in range(max(prefixes)+1)]
    for prefix in prefixes:
        vector[prefix] = 1

    return sum([vector[i]*spectrum[i] for i in range(len(vector))])


def PSM_Search(spectrum, proteome, threshold):

    """ Finds best scoring candidate peptide against spectrum from proteome"""

    candidates = []    
    
    # set bounds for candidate peptide size from min/max AA masses

    low = (len(spectrum) - 1) // max(mass_aa.keys()) + 1
    hi =  (len(spectrum) - 1) // min(mass_aa.keys()) + 1    

    # sliding window samples peptide substrings from proteome

    for i in range(len(proteome)-low+2):
        for j in range(i + low, i+hi+1): 
            peptide = proteome[i:j]
            
            # verify peptide mass matches spectrum mass, add to list

            if len(spectrum)-1 == sum(aa_mass[aa] for aa in peptide):
                if peptide not in candidates:
                    candidates.append(peptide)
    
    # Evauluate candidates, ret best scoring peptide if score> threshold

    best_score = threshold
    match =''
    for peptide in candidates:
        if Score(peptide, spectrum) > best_score:
            best_score = Score(peptide, spectrum)
            match = peptide
    if match:
        return match
#  


# MAIN -----------------
   
with open('data/dataset_11866_5.txt','r') as infile:
   
  # data is lines of >spectra, line of >proteome, line of >threshold

    lines = [line.strip() for line in infile.readlines()]
    threshold = int(lines.pop())
    proteome = lines.pop()

  # **add a root node to all spectra**    

    spectra = [[0]+[int(i) for i in line.split()] for line in lines]
    del(lines)


with open("results/dataset_11866_psm.txt",'w') as file:

  # Call PSM search for each spectrum from input set of spectra

    results = []
    for spectrum in spectra:
        results.append(PSM_Search(spectrum, proteome, threshold))
    for result in results:
        if result:
            file.write(result+'\n')
            print(result)





 


