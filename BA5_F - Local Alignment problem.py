#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 21:00:51 2020

@author: jasonmoggridge

LOCAL ALGINMENT PROBLEM

    * return best scoring alignment of subseq of a,b
    * use PAM250
    * add edge to source for every S(i,j) weight 0 to truncate strings
    * S(sink) is the max score of entire matrix to fly to best seq
    * Truncate alignment based on edges from sink to best to source

"""



def Local_Alignment(a,b):
    
    #Build PAM250 scoring matrix
    
    def PAM250():
        
        with open("/Users/jasonmoggridge/Desktop/PAM250.txt", 'r') as file:
            AA = file.readline().strip().split('  ')
            arrays = []
            for _ in AA:
                array = list(file.readline().strip().split(' '))
                while '' in array:
                    array.remove('')
                array = list(int(i) for i in array[1:])
                arrays.append(array)
        
        pam250 = {}
        for i in range(len(AA)):
            for j in range(len(AA)):
                pam250[(AA[i],AA[j])] = arrays[i][j]
                
        return pam250
    
    
    """ This function walks the manhattan grid and stores the longest path
    to point (i,j) in the Scores matrix """
    
    def Scoring_and_Backtracks(a,b,scoring, sigma):
    
        S = [[0 for _ in range(len(b)+1)]for _ in range(len(a)+1)]  
        S[0][0] = 0

        Pointers = [[0 for _ in range(len(b)+1)]for _ in range(len(a)+1)]    
        Pointers[0][0] = 0
        best = 0; best_node = (0,0)

        for i in range(1, len(a) +1):
            S[i][0] = 0            # penalize all the gaps in first col
            Pointers[i][0] = -1

        for i in range(1, len(b) +1):       # penalize all gaps in first row
            S[0][i] = 0
            Pointers[0][i] = -1

        for i in range(1, len(a) +1):       # iter over rest of matrix   
            for j in range(1, len(b) +1):

                S[i][j] = max(0,                # pick 1/4 moves+source
                 S[i-1][j] - sigma,
                 S[i][j-1] - sigma,
                 S[i-1][j-1] + scoring[(a[i-1], b[j-1])])

                if S[i][j] > best:
                    best = S[i][j]
                    best_node = (i,j)

                if S[i][j] == 0:
                    Pointers[i][j] = -1

                elif S[i][j] == S[i][j-1] - sigma:
                    Pointers[i][j] = 1

                elif S[i][j] == S[i-1][j] - sigma:
                    Pointers[i][j] = 0

                elif S[i][j] == S[i-1][j-1] + scoring[(a[i-1], b[j-1])]:
                    Pointers[i][j] = 2
        
        if S[len(a)][len(b)] < best:
            S[len(a)][len(b)] = best
            Pointers[i][j] = best_node
        
        return S, Pointers
    
    
    """ Function rebuilds the alignment by walking along
    backpointers through the grid, adding a base or gap as indicated by 
    [down: 0, right: 1, diag: 2]    """
    
    def Build_alignment(Pointers, a, b): 
        # walks through pointers in rev topo order to build aligned sequences
    
        if type(Pointers[len(a)][len(b)]) is tuple:
            (i, j) = Pointers[len(a)][len(b)]
        else: 
            i=len(a); j =len(b)
        
        align1 = ''; align2 = ''
    
        while i or j:           
    
            if Pointers[i][j] == -1:          # pointer to (0,0) stop alignment here.
                break
    
            elif Pointers[i][j] == 0:          # pointer say go up a row
                align1 += a[i-1]                 # add the base to align1  
                align2 += '-'                    # add a gap to align2
                i -= 1                           # decrement row number
    
            elif Pointers[i][j] == 1:          # pointer says up move, same as above for seq b
                align1 += '-'
                align2 += b[j-1]
                j -= 1
    
            elif Pointers[i][j] == 2:          # diagonal move, respective bases to each align
                align1 += a[i-1]
                align2 += b[j-1]
                i -= 1; j -= 1
        
        return align1[::-1], align2[::-1]       # reverse each sequence
        #
        
        
    # LOCAL ALIGNMENT Wrapper function   
    Scores, Pointers = Scoring_and_Backtracks(v, w, PAM250(), 5)  
    local_align_score = Scores[-1][-1] 
    return list(Build_alignment(Pointers, v,w)), local_align_score

########



with open("/Users/jasonmoggridge/Desktop/dataset_247_10.txt",'r') as infile:
    v = infile.readline().strip()
    w = infile.readline().strip()

align, local_align_score = Local_Alignment(v,w)

print(local_align_score,'\n',align[0],'\n', align[1])

with open("/Users/jasonmoggridge/Desktop/rosalind_BA5F_out.txt",'w') as outfile:
    outfile.write(str(local_align_score) + '\n')
    outfile.write(align[0] + '\n')
    outfile.write(align[1])
    outfile.close()
