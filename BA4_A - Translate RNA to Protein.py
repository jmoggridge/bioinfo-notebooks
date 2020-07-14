#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 22:11:47 2019
@author: jasonmoggridge
http://rosalind.info/problems/ba4a/


Protein Translation Problem
        Translate an RNA string into an amino acid string.
        
        Given: An RNA string Pattern.
        
        Return: The translation of Pattern into an amino acid string Peptide.
        

"""

def TranslateRNA(rna):

    """ Creates codon_table, then parses RNA into Codons for lookup,
    aminoacids are appended to the string protein and returned
    """
    rna = rna.lower().replace('\n', '').replace(' ', '')
    
    ### codon table ###
    
    bases = [ 'u', 'c', 'a', 'g']
    codons = [a+b+c for a in bases for b in bases for c in bases]
    aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    codon_table = dict(zip(codons, aminoacids))
    
    ### codon lookup ### 
    
    pos = 0
    protein = ''
    while pos < len(rna)-2:
        codon = rna[pos:pos+3]
        for key in codon_table:
            if codon == key:
                if codon_table[key] != '*':
                    protein = protein + codon_table[key]
                    pos +=3
                else:
                    pos +=3
                    break                    
    return (protein)

infile = open('/Users/jasonmoggridge/Desktop/rosalind_ba4a.txt', 'r')
rna = infile.read()

# test case # rna ="AUGCGUA"

prot = TranslateRNA(rna)
# print(str(prot))

output = open("/Users/jasonmoggridge/Desktop/rna_to_prot_out.txt", 'w')
output.write(prot)
output.close()
