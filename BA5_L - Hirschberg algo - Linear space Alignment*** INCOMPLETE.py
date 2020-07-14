#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 21:27:20 2020

@author: jasonmoggridge

Global Alignment in Linear Space Problem:
    
    Find the highest-scoring alignment between two strings using a scoring matrix 
    in linear space.

Given:
    Two long amino acid strings (of length approximately 10,000).

Return: 
    The maximum alignment score of these strings, followed by an 
    alignment achieving this maximum score. Use the BLOSUM62 scoring matrix
    and indel penalty σ = 5.



LinearSpaceAlignment(v, w, top, bottom, left, right)

    # Base cases - zero-area rectangles: return straight indel seq

    if left = right
        output path formed by bottom − top vertical edges
    if top = bottom
        output path formed by right − left horizontal edges

    middle ← ⌊ (left + right)/2⌋
    midEdge ← MiddleEdge(v, w, top, bottom, left, right)
    midNode ← vertical coordinate of the initial node of midEdge 

    LinearSpaceAlignment(v, w, top, midNode, left, middle)

    output midEdge
    if midEdge = "→" or midEdge = "↘"
        middle ← middle + 1
    if midEdge = "↓" or midEdge ="↘"
        midNode ← midNode + 1 

    LinearSpaceAlignment(v, w, midNode, bottom, middle, right)

"""

#def Blosum62():
#    B62 = {}
#    with open("/Users/jasonmoggridge/Dropbox/Rosalind/textbook_track/Course3_Alignments/BLOSUM62.txt", 'r') as file:
#        AA = file.readline().strip().split('  ')
#        arrays = []
#        for _ in AA:
#            array = list(file.readline().strip().split(' '))
#            while '' in array:
#                array.remove('')
#            array = list(int(i) for i in array[1:])
#            arrays.append(array)
#    for i in range(len(AA)):
#        for j in range(len(AA)):
#            B62[(AA[i],AA[j])] = arrays[i][j]
#    return B62
#
#
#def Mid_Edge(v,w):
#
#    def Paths_to_mid(v,w):
#        S = [[0,-sigma]] + list([-sigma*i,-float('inf')] for i in range(1,len(v)+1))    
#        for j in range(1,len(w)+1):
#            for i in range(len(v)+1):
#                match = mu[(v[i-1], w[j-1])]
#                S[i][j%2] = max(\
#                 S[i-1][j%2] - sigma,\
#                 S[i][(j-1)%2] - sigma,\
#                 S[i-1][(j-1)%2] + match)
#        return list(S[i][j%2] for i in range(len(v)+1))
#    
#    mid_j = len(w)//2
#    From_Source = Paths_to_mid(v, w[:mid_j])
#    To_Sink = Paths_to_mid(v, w[len(w):mid_j:-1])
#    best = -float('inf')
#    for i in range(len(v)+1):
#        len_ipath = From_Source[i] + To_Sink[i]
#        if len_ipath > best:
#            best = len_ipath
#            mid_i = i
#
#    if mu[(v[mid_i], w[mid_j])] < -sigma:
#        return  mid_i, (0, 1)
#    else:
#        return mid_i, (1, 1)


def Linear_space_align( v, w, top, bottom, left, right):
    
    print('\n rectangle t/b/l/r:', top,bottom,left,right)
    
    if left == right:
        print('\tleft == right', left, right)
        string =''
        for i in range(top, bottom):
            string += 'D' 
        path.append(string)
        print(string)
        return
    
    if top == bottom:
        print('\tleft == right', left, right)
        
        string = ''
        for j in range(left, right):
            string += 'H'
        path.append(string)
        print(string)
        return        
    
    mid_j = (left+right)//2
    mid_i, edge = Mid_Edge(v[top:bottom],w[left:right])
    
    print('\tNew mid_ij, edge', mid_i, mid_j,edge)
    
    Linear_space_align( v, w, top, mid_i, left, mid_j)
    
    if edge == (0,1):
        path.append('H')
        print('\t\t H')
    elif edge == (1,1):
        path.append('M')
        print('\t\t M')
        
    mid_i += edge[0]
    mid_j += edge[1]
    
    Linear_space_align( v, w, mid_i, bottom, mid_j, right)
    
    ###     
    
path =['$']
with open("/Users/jasonmoggridge/Desktop/rosalind_test.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()

mu = dict(Blosum62())
sigma = 5
Linear_space_align(v,w,0,len(v)+1, 0, len(w)+1)
