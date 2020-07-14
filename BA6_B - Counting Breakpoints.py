"""
BA6_B- 
Number of Breakpoints Problem: 
    Find the number of breakpoints in a permutation.

Input: A permutation.
Output: The number of breakpoints in this permutation.
Code Challenge: Solve the Number of Breakpoints Problem.

Extra Dataset

Sample Input:

+3 +4 +5 -12 -8 -7 -6 +1 +2 +10 +9 -11 +13 +14
Sample Output:

8
"""

with open("/Users/jasonmoggridge/Desktop/rosalind_ba6b.txt",'r') as infile:
    line = infile.readline().strip('(').strip(')')
    P = [int(i) for i in line[:-1].split(' ')]
    infile.close()

def breakpoints(P):
    P = [0] + P + [len(P)]    
    breakpoints = 0
    for k in range(len(P)-1):
        if P[k+1]- P[k] != 1:
            breakpoints +=1
    print('breakpoints:', breakpoints)
    return breakpoints
