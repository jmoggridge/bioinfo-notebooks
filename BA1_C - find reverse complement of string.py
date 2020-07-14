#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 19:13:56 2019
@author: jasonmoggridge

Find the Reverse Complement of a String 
←→:
    In DNA strings, symbols 'A' and 'T' are complements of each other,
    as are 'C' and 'G'. Given a nucleotide p, we denote its complementary
    nucleotide as p. The reverse complement of a DNA string (Pattern = p1…pn)
    is the string Pattern = pn … p1 formed by taking the complement of each
    nucleotide in Pattern, then reversing the resulting string.

For example, the reverse complement of Pattern = "GTCA" is Pattern = "TGAC".

Reverse Complement Problem

Find the reverse complement of a DNA string.

Given: A DNA string Pattern.

Return: Pattern, the reverse complement of Pattern.

Sample Dataset
AAAACCCGGT

Sample Output
ACCGGGTTTT


"""


# Reads seqs -> and builds reverse complement rc <-backwards from     
def ReverseComp(read):
    
    complement = {'A':'T','C':'G','G':'C','T':'A'}
    read_rc = ''
    for base in read:
        read_rc = complement[base] + read_rc
    
    return read_rc
    ###
    
print("\nInput your DNA string:\n\t")
DNA = input()
RC = ReverseComp(DNA)

print('\n\nReverse Complement of input DNA:\n\n', RC)

with open("/Users/jasonmoggridge/Desktop/rosalind_ba1c_out.txt",'w') as op:
    op.write(RC)
    op.close()