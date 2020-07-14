#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 03:49:41 2020

@author: jasonmoggridge
"""

def Mass_Table():  
    mass_table ={}
    with open('integer_mass_table.txt', 'r') as file:
        for line in  file.readlines():
            [aa,mass] = line.strip().split()
            mass_table[int(mass)] = str(aa)
    return mass_table
mass_table = Mass_Table()
##

def Graph(S):
    graph = {S[i]:[] for i in range(len(S))}
    for i in range(len(S)-1):
    
        for j in range(i, len(S)):
     
            delta_ij = S[j] - S[i]
            if delta_ij in mass_table.keys():
                
                a_ij = mass_table[delta_ij]
                graph[S[i]].append((S[j], a_ij))
    return graph



# Main

with open('rosalind_ba11a.txt','r') as infile:
    S = [0]+[int(i) for i in infile.readline().strip().split()]

graph = Graph(S)

with open('results/spectrum_graph_out.txt','w') as outfile:
    string =''
    for mass_i in S:
        for j in graph[mass_i]:
            string += str( str(mass_i) + '->' + str(j[0]) + ':' + str(j[1])+'\n')
    string.strip("\n")
    outfile.write(string)