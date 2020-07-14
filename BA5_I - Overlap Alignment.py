#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 00:08:56 2020

@author: jasonmoggridge

Code Challenge: Solve the Overlap Alignment Problem.

Input:
    Two strings v and w, each of length at most 1000.

Output:
    The score of an optimal overlap alignment of v and w,
 followed by an alignment of a suffix v' of v and a prefix w'
 of w achieving this maximum score. Use an alignment score in
 which matches count +1 and both the mismatch and indel
 penalties are 2.


"""


def Overlap_Alignment(v,w):

    #  Use an alignment score in which matches count +1
    # and both the mismatch and indel penalties are 2.

    sigma = 2
    mu = sigma

    # Setup a pointers matrix for backtracking moves.

    Pointers = [[False for _ in range(len(w)+1)]for _ in range(len(v)+1)]
    Pointers[0][0] = 0

    # Give free ride through prefix of v in first column of S matrix
    # but penalize any initial gaps in w in first row

    S = [[0 for _ in range(len(w)+1)]for _ in range(len(v)+1)]
    S[0][0] = 0
    for i in range(1, len(v) +1):
        S[i][0] = 0
        Pointers[i][0] = 0
    for i in range(1, len(w) +1):
        S[0][i] = -sigma*i
        Pointers[0][i] = 1

    # Fill out scores for entire matrix(v,w)

    for i in range(1, len(v) +1):
        for j in range(1, len(w) +1):

            match = -mu
            if v[i-1] == w[j-1]:
                match = 1
            S[i][j] = max(
                     S[i -1][j] - sigma,
                     S[i][j-1] - sigma,
                     S[i-1][j-1] + match)

            # pointer code: down 0, right 1, diag 2

            if S[i][j] == S[i-1][j] - sigma:
                Pointers[i][j] = 0
            elif S[i][j] == S[i][j-1] - sigma:
                Pointers[i][j] = 1
            elif S[i][j] == S[i-1][j-1] - mu:
                Pointers[i][j] = 2
            elif S[i][j] == S[i-1][j-1] + 1:
                Pointers[i][j] = 2

    # Figure out the best scoring overlap by scanning max(bottom row of S)

    score = -float('inf')
    end = 0

    for j in range(len(w),0,-1):
        if S[-1][j] > score:
            score = S[-1][j]
            end = j

    return Pointers, end, score, S


# walk through pointers in rev topo order from 'end'+ build alignmt

def Build_alignment(Pointers, v, w, end):

    align1 = ''; align2 = ''
    i = len(v)
    j = end

    # Backtrack + Align until seq w runs out, truncate v here.
    while j > 0:

        if Pointers[i][j] == 0: # down
            align1 += v[i-1]
            align2 += '-'
            i -= 1

        elif Pointers[i][j] == 1: # right
            align1 += '-'
            align2 += w[j-1]
            j -= 1

        elif Pointers[i][j] == 2: # diagonal
            align1 += v[i-1]
            align2 += w[j-1]
            i -= 1; j -= 1

    return align1[::-1], align2[::-1]


###


with open("/Users/jasonmoggridge/Desktop/rosalind_ba5i.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()
    infile.close()

Pointers, end, score, S = Overlap_Alignment(v, w)
align = Build_alignment(Pointers, v, w, end)
print(score,'\n',align[0],'\n', align[1])

with open("/Users/jasonmoggridge/Desktop/rosalind_ba5i_out.txt",'w') as outfile:
    outfile.write(str(score) + '\n')
    outfile.write(align[0] + '\n')
    outfile.write(align[1])
    outfile.close()
