#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 23:44:27 2019

@author: jasonmoggridge
     
    """

def genes_in_seq(DNA,peptide):
    
    
    # Builds codon table
    def codontable():
        bases = [ 'T', 'C', 'A', 'G']
        codons = [a+b+c for a in bases for b in bases for c in bases]
        aminoacids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
        codon_table = dict(zip(codons, aminoacids))
        codon_table.update(zip(aminoacids,codons))
        return codon_table
        
    # Reads seqs -> and builds reverse complement rc <-backwards from     
    def ReverseComp(read):
        complement = {'A':'T','C':'G','G':'C','T':'A'}
        read_rc = ''
        for base in read:
            read_rc = complement[base] + read_rc
        return read_rc
    ##
        
    def TranslateDNA(dna):
        i = 0
        protein = ''
        while i < len(dna)-2:
            codon = dna[i:i+3]
            protein += codons[codon]
            i+=3
        return protein    
    ##
    
    def find_genes(peptide, DNA):
    
        sequences = []
        for shift in range(3):
            protein = TranslateDNA(DNA[shift:])
            for i in range(len(protein) - len(peptide)+1):
                if protein[i:i+len(peptide)] == peptide:
                    sequences.append(DNA[ 3*i + shift : 3*i + 3*len(peptide) + shift])
    
        return sequences
    ##

    codons = codontable() 
    sequences = find_genes(peptide, DNA)
    sequences += list(ReverseComp(seq) for seq in find_genes(peptide, ReverseComp(DNA)))

    return sequences
### 
    
file = open("/Users/jasonmoggridge/Desktop/rosalind_ba4b.txt", 'r')
DNA = str(file.readline().strip())
peptide = str(file.readline().strip())
###

sequences = genes_in_seq(DNA, peptide)
wr = open("/Users/jasonmoggridge/Desktop/rosalind_test_output.txt", 'w')

for seq in sequences:
    wr.write(seq+'\n')
wr.close()
    

