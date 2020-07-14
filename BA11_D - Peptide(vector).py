#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 07:17:25 2020

@author: jasonmoggridge
"""

def Mass_Table():  

    mass_aa = {}
    aa_mass = {}

    with open('integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_aa[int(mass)] = str(aa)
            aa_mass[str(aa)] = int(mass)

    return mass_aa, aa_mass
mass_aa, aa_mass = Mass_Table()
#
def Peptide(vector):
    
    peptide = ''
    sum_mass = 0
    for mass in range(len(vector)):
        if vector[mass]:
            delta = mass + 1 - sum_mass
            peptide += mass_aa[delta]
            sum_mass += delta
    return peptide


#
#def Vector(peptide):
#    prefixes = [sum(aa_mass[aa] for aa in peptide[:i]) for i in range(1,len(peptide)+1)]
#    vector = [0 for _ in range(max(prefixes))]
#    for prefix in prefixes:
#        vector[prefix-1] = 1
#    return vector
#
# #vector = Vector(peptide)



# test data

#with open('rosalind_test.txt', 'r') as infile:
   #mass_aa = {4:'X', 5:'Z'}
   #aa_mass = {'X':4, 'Z':5}

#with open('peptide_vector_to_peptide.txt', 'r') as infile:
#    infile.readline()
#    vector = [int(i) for i in infile.readline().strip().split()]
#    infile.readline()
#    solution = infile.readline().strip()
#peptide = Peptide(vector)
#print(peptide)
#print('correct?:', peptide == solution)

### Graded data
with open('dataset_11813_8.txt', 'r') as infile:
    vector = [int(i) for i in infile.readline().strip().split()]
peptide = Peptide(vector)
print(peptide)

