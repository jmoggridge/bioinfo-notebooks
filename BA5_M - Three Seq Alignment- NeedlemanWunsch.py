#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 00:16:40 2020

@author: jasonmoggridge

"""

def Multiple_Sequence_Alignment(v,w,u):
        
    # Scoring: 1 if all bases match, else 0        
    def Score(a,b,c):
        if a == b == c:
            return 1
        return 0

    # this is an array of vectors corresponding to incoming edges
    # for storing backpointers as an integer and for manipulating graph co-ordinates
    # when creating the alignment later on
    moves = [(-1,0,0), (0,-1,0), (0,0,-1), (-1,-1,0), (-1,0,-1), (0,-1,-1), (-1,-1,-1)]    
    
    # S matrix holds dynamic programming info for optimal path length to node_ijk
    # where S_ijk = len(longest path to node_ijk)
    import numpy as np
    S = np.zeros( (len(v), len(w), len(u)), dtype=np.int16)
    
    # P matrix stores the index of vector of the previous move to node_ijk
    P = -1* np.ones( (len(v),len(w),len(u)), dtype=np.int16)
    
    # The pointers are initialized such that indels at the ends are minimized
    # by backpointing vectors to the source being 2 dimension whenever possible
 
    # init all (0,0,k) pointers to (0,0,-1) ie. moves[2]    
    for k in range(1, len(u)): 
        P[0][0][k] = 2
     
    # init all (0,j,0) pointers to (0,-1,0) ie. moves[1]   
    # init all (0, [1:], [1:]) plane to (0,-1,-1) ie. moves[5]   
    
    for j in range(1, len(w)):
        P[0][j][0] = 1  
        for k in range(1, len(u)):
            P[0][j][k] = 5
    

    # for matrix i; set the pointer for (i,0,0) to (-1,0,0) -> moves[0]
    for i in range(1, len(v)):
        P[i][0][0] = 0 
        
        # for row j; set the first pointer in each row to (-1,-1,0): move[3]
        for j in range(1, len(w)):
            P[i][j][0] = 3 
        
            # for column k; first pointer in each col to (-1,0,-1): move[4]
            for k in range(1, len(u)):
                P[i][0][k] = 4
                
                # Store the length of longest path to i,j,k
                S[i][j][k] = max(
                        S[i-1][j-1][k-1] + Score(v[i], w[j], u[k]),
                        S[i-1][j-1][k],
                        S[i-1][j][k],
                        S[i][j-1][k],
                        S[i-1][j][k-1],
                        S[i][j-1][k-1],
                        S[i][j][k-1])
                
                # Retroactively figure out which of 7 moves was the optimal edge
                for move in moves:
                    
                    # co-ords and score of parent node:
                    pi,pj,pk = i + move[0], j + move[1], k + move[2]
                    score = S[pi][pj][pk]
                    
                    # need to account for mis/match scoring for this edge
                    if move == (-1,-1,-1):
                        score += Score(v[i], w[j], u[k])
                    
                    # if this parent + edge is what we stored in Sijk,
                    # that is an optimal parent, store the vector in P_ijk
                    if S[i][j][k] == score:
                        P[i][j][k] = moves.index(move)


    # Now to retrace backpointers to construct aligned sequences
    
    Seqs = [v[1:],w[1:],u[1:]]              
    Aligns = ['','','']
    
    # Global alignment so start backtracking from end of all seqs
    ijk = [len(v)-1, len(w)-1, len(u)-1]    
    
    # continue building alignment until all sequences finished (global align)
    while sum(ijk) > 0: 
        
        # Get backtrace vector
        move = moves[P[ijk[0]][ijk[1]][ijk[2]]]     
        
        # Add a character to each aligned sequence
        # either a gap or base depending whether edge is 0 or -1 in that dimension
        for seq in range(len(Seqs)):
            if move[seq] == 0:                      
                Aligns[seq] += '-'
            elif move[seq] == -1:
                Aligns[seq] += Seqs[seq][ijk[seq]-1]
                ijk[seq] -=1
    
    #print score       
    print(S[-1][-1][-1])
    
    # reverse aligned seq and print
    Aligns = [align[::-1] for align in Aligns]
    for align in Aligns:
        print(str(align))
    return int(S[-1][-1][-1]), Aligns
###

with open("/Users/jasonmoggridge/Desktop/rosalind_ba5m.txt",'r') as infile:
    v = '-'+ infile.readline().strip()
    w = '-'+ infile.readline().strip()
    u = '-'+ infile.readline().strip()                
score, aligns = Multiple_Sequence_Alignment(v,w,u)

        

