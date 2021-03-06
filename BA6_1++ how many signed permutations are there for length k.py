#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 23:52:19 2020

@author: jasonmoggridge
"""

""" How many unique signed permutations of length n are there???"""

def count_signed_permutations(n):

    N = list(2*i for i in range(1,n+1))
    perms =1
    for i in N:
        perms = perms*i
    return perms


n=7
count_signed_permutations(n)





"""Exercise Break: How many different reversals can 
be applied to a permutation of length 100???"""

def how_many_reversals(k):
    return int((k+1)*k/2)

k=100
how_many_reversals(k)




"""How ordered is a permutation???
    Reverssal distance -> min #reversals P -> Q
        *very complicated algorithm*

    GreedySorting(P)
        approxReversalDistance ← 0
        for k = 1 to |P|
            if element k is not sorted
                apply the k-sorting reversal to P
                approxReversalDistance ← approxReversalDistance + 1
            if k-th element of P is −k
                apply the k-sorting reversal to P
                approxReversalDistance ← approxReversalDistance + 1
        return approxReversalDistance                           """
        
def Greedy_sorting(P):
    
    def k_sorting_reversal(P,k):
        
        sign = 1
        if -k in P: sign = -1        
        rev = []
        for i in range(k-1, P.index(k*sign)+1):
            rev = [P[i]*-1] + rev
        return P[:k-1] + rev + P[P.index(k*sign)+1:]
    #
    
    approx_distance = 0
    for k in range(1,len(P)+1):

        if P[k-1] != abs(k):
            P = k_sorting_reversal(P,k)
            approx_distance += 1
            line = ' '.join(('+' if i > 0 else '') + str(i) for i in P)
            print(line)
            outfile.write(line)

        if P[k-1] != k:
            P[k-1] = k
            approx_distance += 1
            line = ' '.join(('+' if i > 0 else '') + str(i) for i in P)
            print(line)
            outfile.write(line)
           
    return approx_distance
##

outfile = open("/Users/jasonmoggridge/Desktop/greedy_sorted.txt",'w')
with open("/Users/jasonmoggridge/Desktop/dataset_286_4.txt",'r') as infile:
    P = [int(i) for i in infile.readline().strip(')').strip('(').split(' ')]
print('\n')
approx_d = Greedy_sorting(P)


outfile.close()