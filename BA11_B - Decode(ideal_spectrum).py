#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 05:24:54 2020

@author: jasonmoggridge


#
#DecodingIdealSpectrum(Spectrum)
#     construct Graph(Spectrum)
#     for each path Path from source to sink in Graph(Spectrum)
#          Peptide ← the amino acid string spelled by the edge labels of Path
#          if IdealSpectrum(Peptide) = Spectrum
#                return Peptide

"""

def Decode(Spectrum):
    
    # decode(spectrum) subroutines: Mass_aa dicts; IdealSpectrum; dfs_peptides; decode_wrapper

    def Mass_Table():  
    
        mass_aa = {}
        aa_mass = {}
    
        with open('integer_mass_table.txt', 'r') as file:
            for line in  file.readlines():
                [aa,mass] = line.strip().split()
                mass_aa[int(mass)] = str(aa)
                aa_mass[str(aa)] = int(mass)
    
        return mass_aa, aa_mass
    
    #
    
    def IdealSpectrum(peptide):
    
        spectrum = []
        for i in range(len(peptide)):
            spectrum.append(sum(aa_mass[aa] for aa in peptide[:i]))
            spectrum.append(sum(aa_mass[aa] for aa in peptide[i:]))
        return sorted(spectrum)
        #
        
    def Graph(S):
        tree = {S[i]:[] for i in range(len(S))}
        for i in range(len(S)-1):
        
            for j in range(i, len(S)):
         
                delta_ij = S[j] - S[i]
                if delta_ij in mass_aa.keys():
                    
                    a_ij = mass_aa[delta_ij]
                    tree[S[i]].append((S[j], a_ij))
        return tree
        #
    
    def dfs_peptides(graph):
        
        def dfs_search(graph, v, path, explored):
           
            explored[v] = True
            for w in graph[v]:
                if w[0] == sink:
                    if str(path+w[1])[::-1] not in paths:
                        paths.append(path + w[1])
                else:
                    dfs_search(graph, w[0], path+w[1], explored)
    
        sink = max(graph.keys())
        paths = []
        explored = {key:False for key in graph.keys()}
        dfs_search(graph, 0, '', explored)
        return paths
        #

    # Decode(spectrum) function:
    mass_aa, aa_mass = Mass_Table()
    graph = Graph(Spectrum)
    peptides = dfs_peptides(graph)
    for peptide in peptides:
        if IdealSpectrum(peptide) == Spectrum:
            return peptide
#


#Main

with open('rosalind_ba11b.txt','r') as infile:
    Spectrum = [0]+[int(i) for i in infile.readline().strip().split()]

peptide = Decode(Spectrum)
print('\n\nIdeal spectra match for peptide (or its reverse):\n\n\t', peptide)
