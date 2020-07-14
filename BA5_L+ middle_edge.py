#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 01:08:10 2020

@author: jasonmoggridge
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 16:17:52 2020

@author: jasonmoggridge

Find a middle edge in the alignment graph in linear space.

Given: Two amino acid strings.

Return: A middle edge in the alignment graph of these strings,
 where the optimal path is defined by the BLOSUM62 scoring matrix 
 and a linear indel penalty equal to 5. Return the middle edge in the 
 form “(i, j) (k, l)”, where (i, j) connects to (k, l).
 

A space-efficient seq algorithm

    - use two column structure to store (_,j)(_,j-1) values for longest path
        - split into two subproblems for tracing best alignment
    - identify middle node
    - divide and conquer

Code Challenge: Solve the Middle Edge in Linear Space Problem (for protein strings).

Input: Two amino acid strings.

Output: A middle edge in the alignment graph in the form "(i, j) (k, l)",
 where (i, j) connects to (k, l). To compute scores, use the BLOSUM62
 scoring matrix and a (linear) indel penalty equal to 5.
 
 Middle Edge in Linear Space Problem

 
 """



def Blosum62():
    B62 = {}
    with open("/Users/jasonmoggridge/Dropbox/Rosalind/textbook_track/Course3_Alignments/BLOSUM62.txt", 'r') as file:
        AA = file.readline().strip().split('  ')
        arrays = []
        for _ in AA:
            array = list(file.readline().strip().split(' '))
            while '' in array:
                array.remove('')
            array = list(int(i) for i in array[1:])
            arrays.append(array)
    for i in range(len(AA)):
        for j in range(len(AA)):
            B62[(AA[i],AA[j])] = arrays[i][j]
    return B62



def iPaths(v,w):

    S = [[0,-sigma]] + list([-sigma*i,-float('inf')] for i in range(1,len(v)+1))    
    for j in range(1, len(w) + 1):
        for i in range(1, len(v) + 1):
#
#            match = mu[(v[i-1], w[j-1])]
            match = -mu
            if v[i-1] == w[j-1]:
                match = mu
            S[i][j%2] = max(\
             S[i-1][j%2] - sigma,\
             S[i][(j-1)%2] - sigma,\
             S[i-1][(j-1)%2] + match)

    return list(S[i][j%2] for i in range(len(v)+1))



def Mid_node( v, w, top, bottom, left, right):
    
    print('\n\nMidNode - v, w, top, bottom, left, right')
    print(v,w,top,bottom,left, right)
    
    mid_j = (left + right) //2
    print('Mid_j:', mid_j)
    
    v = v[top:bottom]
    w = w[left:right]
    From_Source = iPaths( v, w[: mid_j])
    
    print('\nSlice of v: (top:bottom):', v, top, bottom)
    print('Slice of w: (left:mid_j+1):', w[:mid_j], left, mid_j)
    print('From_source:', From_Source)
    
    rev_v = v[::-1]
    print('\n\nReversed v for To_Sink:', rev_v)
    rev_w = w[::-1]
    print('Reversed w for To_Sink:', rev_w)
    
    To_Sink = iPaths( rev_v[top:bottom], w[:mid_j+1])
    To_Sink.reverse()
    print('\nSlice of v rev: (bottom: top):', str(rev_v[top: bottom]), top, bottom)
    print('Slice of w rev: right -> mid_j+1', str(rev_w[: mid_j+1]), right, mid_j+1)
    print('To_sink:', To_Sink)

    best = -float('inf')
    print('\n\tscoring: ', best)
    for i in range(top, bottom):
        print('\t\tI:', i)
        len_ipath = From_Source[i] + To_Sink[i]
        print('\t\t\tlen_ipath = ', len_ipath, '= From_Sr:', From_Source[i],'+ to sink:', To_Sink[i])
        if len_ipath >= best:
            print('\t\tscoring: ', best, '@ I:' , i, 'len:',len_ipath)
            best = len_ipath
            mid_i = i + top

#    score = max(-sigma, mu[(v[mid_i], w[mid_j])])

    return  (mid_i, mid_j)


#with open("/Users/jasonmoggridge/Desktop/rosalind_test.txt",'r') as infile:
#    v = infile.readline().strip()
#    w = infile.readline().strip()
##mu = dict(Blosum62())


v = 'TTTACC'
w = 'ATCT'
#
v = 'ACTTAATT'
w = 'T'

n = len(v)
m =len(w)
mu = 1
sigma = 2

node =  Mid_node( v, w, 0, len(v)+1, 0, len(w)+1)

