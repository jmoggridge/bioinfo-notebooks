#!usr/bin/env python

""" ManberMyers suffix array constructor"""


import numpy as np


### Subroutines

def inv_arr(SA):
    """Invert suffix array for quick lookups"""
    ISA = list(range(len(SA)))
    for i in range(len(SA)):
        ISA[SA[i][0]] = i 
    return ISA

def update_ranks(SA,n):
    """Make entries with matching prefixes have the same rank after the previous sorting step"""
    # stores last entry from ranks(i-1,i-1+h)
    prev = list(SA[0][1:])  
    # rank of lex.first suffix is zero
    SA[0][1] = 0
    
    for i in range(1,n):
        if np.array_equal(SA[i][1:], list(prev)): 
            # assign same rank as i-1 if same h-prefix 
            SA[i][1] = SA[i-1][1]
        else:
            # rank(i,i+h)_i is > prev(i-1, i-1+h), tf. rank[i] = rank[i-1]+1
            prev = list(SA[i][1:])
            SA[i][1] = SA[i-1][1]+1     
    return SA

def pair_ranks(SA, ISA, h, n): # might need to add counts here...? don't understand,,
    """find matching key for each S(i) -> S(i+h)"""
    for i in range(n):
        j = SA[i][0] + h
        if j >= n:
            SA[i][2] = -1
        else:
            SA[i][2] = SA[ISA[j]][1]
    return SA

def two_radix(SA):
    """radix sort of suffix array (suffix, rank_i, rank_i+h) matrix"""
    for column in range(2,0,-1):
        SA = SA[np.argsort(SA[:,column])]
    return SA

def initial_SA(text,n):
    """Initialize suffix matrix SA[ [0:suffix ranked i],[1:rank of suffix],[2:rank of (i,i_h)-paired suffix]"""
    SA =[[i,rank,0] for i, rank in zip(range(n), (ord(i) for i in text))]
    SA = np.array(SA)
    SA = SA[np.argsort(np.array(SA)[:,1])]
    return SA

def BWT(SA,text):
    """SA_i starts with i-th character in text -> bwt is the (i-1)th character in text"""
    return ''.join(text[row[0]-1] for row in SA)
        
### Wrapper

def manber_myers(text):
    
    if text[-1] != '$':
        text += '$'
    
    h=1
    n = len(text)
    SA = initial_SA(text,n)  
    ISA = inv_arr(SA)
    
    while h < 2*n+1: 
        SA = update_ranks(SA,n)
        SA = pair_ranks(SA,ISA,h,n)
        SA = two_radix(SA)    
        h=2*h
        ISA = inv_arr(SA)
    bwt = BWT(SA, text)
    SA = list(int(SA[i][0]) for i in range(n))
    return SA, bwt


################################################################################################################################

### Table display for Ipython notebooks
    
#def SA_table(SA,ISA,text):
#    """tabulate function for Manber-Myers algorithm"""
#    sa = np.ndarray.tolist(SA)
#    for i in range(len(sa)):
#        si = text[sa[i][0]:] + text[:sa[i][0]]
#        
#        sa[i]= [si, sa[i][0], ISA[i]] + sa[i][1:]
#    head = ['i','suffix','S[i]','ISA[i]','Rank[i]', 'Rank[i-2h]']
#    display(HTML(tabulate.tabulate(sa,headers=head, tablefmt='html', showindex="always")))  #
#    return



## Test unit classix:
    
text = 'panamabananas$'
SA,bwt = manber_myers(text) 
print(text)
print(bwt)
print(', '.join(str(s) for s in SA))




text = 'AACGAATAAGAACGAGATCGGTAGA$'
SA,bwt = manber_myers(text) 
print(text)
print(bwt)
print(', '.join(str(s) for s in SA))


with open("/Users/jasonmoggridge/Desktop/rosalind_ba9g_1_dataset.txt") as infile:
    text = infile.readline().strip()


with open("/Users/jasonmoggridge/Desktop/rosalind_ba9g_1_output.txt") as infile:
    solution = infile.readline()
    sol_list = solution.split(', ')
###

n=len(text)
SA, bwt = MM.manber_myers(text)
SA = list(str(s) for s in SA)


with open('/Users/jasonmoggridge/Desktop/rosalind_out.txt', 'w') as outf:
    outf.write(', '.join(str(s) for s in SA))
    outf.write('\n')
    
print(SA == sol_list)

for i in SA[:100]:
    print(text[int(i):int(i)+10])