#!usr/bin/env python

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import itertools
import time
import operator
import sys

class BurrowsWheeler: 
    """BWT & FM-index object for a text sequence"""

    def __init__(self, seq, **k): # find ideal k and alphabet
        """BW(seq,k) generates FM-index:\n BW.bwt, FirstOccurence, Checkpoint array, Partial suffix array \
        k- is the factor by which the FM-index (expt. bwt) is reduced from it's original size"""
        
        if seq [-1] != "$":
            seq += "$"
        self.seq = seq
        
        self.alpha = sorted(set(seq))
        if not k:
            k = 5
            print('setting default k')
        self.k = k
        
        #original naive
        start = time.time()
        self.bwt1, self.psa1 = self.bw_transform()   # Burrows-Wheeler transform and k-suffix array
        print('BWT naive sort: time (s) =',time.time()-start)
        
#         #radix
#         start = time.time()
#         self.bwt, self.psa = self.bwt_radix()   # new sort to try
#         print('BWT radix-like: time (s) =',time.time()-start)
        
        #manber myers
        start = time.time()
        self.bwt, self.psa = self.manber_myers()
        print('ManberMyers: time (s) =',time.time()-start)

        
        self.first_occur = self.first_occurence()  # FirstOcc[symbol] pos in text
        self.checkpt = self.get_checkpoint()       # Checkpoint array


    def bwt_radix(self):
        #shouldn't be using numpy whereever possibly...
        
        """radix-sort-like block encoding bwt from text, yields bwt, partial suffix array"""
        
        suffix_array = np.array(range(len(self.seq)), dtype=np.int)
        
        for x in range(len(self.seq)-1,-1,-1):
            rotation = self.seq[x:]+ self.seq[:x]  
            new = [[] for _ in self.alpha]
        
            for i, suffix in enumerate(suffix_array):
                symbol = rotation[suffix_array[i]]
                new[self.alpha.index(symbol)].append(suffix)
            suffix_array = np.array([item for sublist in new for item in sublist], dtype=np.int)

        psa = [-1] * ((len(self.seq)-1)//k + 1)
        bwt = ''.join(self.seq[s-1] for s in suffix_array)
        for suf in suffix_array:
            if suf % k == 0:
                psa[suf//k] = np.where(suffix_array == suf)[0][0]

        return bwt, psa
        
    
    def manber_myers(self):
        start = time.perf_counter()
        n = len(self.seq)
        # Pos array sorted by first symbol in (Pos[i], Prm[i]) = i, ascii-score(text[i])
        Pos = [[i, ord(symbol)] for i, symbol in enumerate(self.seq)]
        Pos.sort(key = operator.itemgetter(1))
        (Pos, Prm) = (list(i) for i in list(zip(*Pos)))

        # init boolean arr
        Bh, B2h = ([False for _ in range(n)] for _ in range(2))
        Bh.append(True)

        # set degenerate ranks for suffixes with same first chars
        temp = Prm[0]       
        Prm[0] = 0
        Bh[0] = True        
        for i in range(1, n):
            if Prm[i] == temp:
                Prm[i] = Prm[i-1]
            else:
                temp = Prm[i]
                Prm[i] = i
                Bh[i] = True

        #init count array
        Count = [0 for _ in range(n)]
        # init setup done


        # while loop
        h = 1

        while h<n:
            # update Prm <- inverse(Pos)
            for i in range(n):
                Prm[Pos[i]] = i

            #reset count array
            Count = [0 for _ in range(n)]

            # set Prm[d] <- Prm[head of bucket[d]]
            head = 0
            for i in range(1, n):
                if Bh[i]:
                    head = Prm[Pos[i]]
                else:
                    Prm[Pos[i]] = head   

            ##  update the suffix that ends this round... ? 
            # not sure if this is actually? doing? anything?
            d = n-h
            e = Prm[d]
            Prm[d] = e + Count[e]
            B2h[Prm[d]] = True

            # For each '=h bucket', update paired keys, set Bh for new heads that arise...
            left = 0
            c = 0
            while c < n:
                d = Pos[c] - h
                if d >= 0:
                    e = Prm[d]
                    Prm[d] = e + Count[e]
                    Count[e] += 1
                    B2h[Prm[d]] = True
                c += 1
                # Bh[c] ->if True, c is end of h-bucket; fix B2h
                if Bh[c]: 
                    for i in range(left,c):
                        d = Pos[i]-h
                        # B2h(Prm[d]) is a new head, 
                        if d >= 0 and B2h[Prm[d]]: 
                            j = Prm[d] + 1
                            # set others in bucket below d to false, end at next head (Bh[j]=Tr)
                            while True:
                                if Bh[j] or not B2h[j]:
                                    break
                                B2h[j] = False
                                j += 1
                # set the left border for the new bucket
                    left = c

            # Pos vector <- updated from Prm
            for i in range(n):
                Pos[Prm[i]] = i
            # Bh <- B2h
            counter=0
            for i in range(n):
                if B2h[i]:
                    Bh[i] = True
                if Bh[i]:
                    counter += 1
            if counter == n:
                print('\nFinished early, at stage h=',h, 'time = ', time.perf_counter()-start)
                break

            h = 2*h

        print('done at h=',h)   

        psa = [-1] * ((len(self.seq)-1)//k + 1)
        bwt = ''.join(self.seq[s-1] for s in Pos)
        Pos = np.array(Pos)
        for suf in Pos:
            if suf % k == 0:
                psa[suf//k] = np.where(Pos == suf)[0][0]

        return bwt, psa
    
    
    def bw_transform(self): 
        """naive bwt with python sorted() on |text|^2 array"""
        cycles = sorted([(self.seq[i:] + self.seq[:i], i) for i in range(len(self.seq))], key= lambda x: x[0])
        bwt = ''.join(cycle[0][-1] for cycle in cycles)
        sa = [len(self.seq) - cycle[0].index('$')-1 for cycle in cycles]
        psa = [-1] * ((len(self.seq)-1)//self.k + 1)
        for suf in sa:
            if suf % self.k == 0:
                psa[suf//self.k] = sa.index(suf)
        return bwt, psa

    def first_occurence(self):
        """dict[symbol] = first index of symbol in sorted(bwt)"""
        first_occur = {}
        for i, symbol in enumerate(sorted(self.bwt)):
            if symbol not in first_occur.keys():
                first_occur[symbol] = i
        return first_occur

    def get_checkpoint(self):
        """Checkpoint matrix for pattern matching, every k-th row of count matrix"""
        alpha = sorted(set(self.seq))
        count = {symbol:[0] for symbol in alpha}
        for i in range(self.k, len(self.bwt), self.k): #self.bwt
            c = Counter(self.bwt[i - self.k:i])
            for s in alpha:
                count[s].append(count[s][-1] + c[s])
        return count   
    
#     # not necessary here...
#     def decode(self):
#         """decodes the BWT to the original string text"""
#         last = [(symbol, self.count[symbol][index]) for index, symbol in enumerate(self.bwt)]
#         first = sorted(last, key = lambda x: x[0])
#         decoded = ''
#         symbol = ('$',0)
#         while len(decoded)<len(last):
#             symbol = first[last.index(symbol)]
#             decoded += symbol[0]
#         return decoded

    def pattern_match(self, pattern):
        
        """BW matching using first occur and count data, constant-time indexing"""
        top = 0
        bottom = len(self.bwt)
        while bottom > top:
            # case- match(es): return index pointers
            if not pattern:
                return (top, bottom)  
            symbol = pattern[-1]
            # case- no match, pattern has symbol not in text
            if symbol not in self.checkpt.keys():
                return None           
            
            # index top & bottom pointers using first_occ + checkpoint + counter
            c_top = self.first_occur[symbol] + self.checkpt[symbol][top//self.k]
            if top % self.k != 0:
                c_top += Counter(self.bwt[top - top%self.k: top])[symbol]
            
            c_bottom = self.first_occur[symbol] + self.checkpt[symbol][bottom//self.k]
            if bottom % self.k != 0:
                 c_bottom += Counter(self.bwt[bottom - bottom % self.k: bottom])[symbol]
            
            # set new pointers and remove symbol from pattern
            (top, bottom) = (c_top, c_bottom)
            pattern = pattern[:-1]
        
        # case- no match, bottom=top  
        return None 
    
    
    def multiple_pattern_match(self, patterns):
        """Returns a set of start indexes of all matches in text (all in one list right now)"""
        
        pattern_windows = list(map(self.pattern_match, patterns)) # PM for each pattern
        
        # recreate First/Last cols -> [(symbol, rank in last)]
        last = []
        for i, s in enumerate(self.bwt):
            if s == '$':
                last.append((s,0))
            else:
                last.append((s, self.checkpt[s][i//self.k]+ Counter(self.bwt[i - i%self.k: i])[s]))
        first = [(self.bwt[i],last[i][1]) for i in range(len(self.bwt))]
        first = sorted(first, key = lambda x: x[0])        
        results = []
        
        # Convert top-bottom windows -> start indexes
        for window in pattern_windows:
            
            # case- no-match-for-pattern 
            if not window:
                results.append(None)
            
            # case- pattern has >=1 match in text
              # window is the (top, bottom) indices in First for each pattern 
              # walk back in text by first/last property until find entry in psa
              # i's start pos. is psa[j] + steps(from i to j)
                
            else:
                starts = []
                for i in range(window[0], window[1]):
                    steps = 0
                    j = i       
                    while True:
                        if j in self.psa: 
                            starts.append(self.psa.index(j)*self.k + steps)
                            break
                        steps+=1
                        j = first.index(last[j]) 
                        
                results.append(sorted(starts))
        return results # list of lists with either [none] or [start positions]


# ManberMyers - algorithm


# def manber_myers(text):
#     start = time.perf_counter()
#     n = len(text)
#     # Pos array sorted by first symbol in (Pos[i], Prm[i]) = i, ascii-score(text[i])
#     Pos = [[i, ord(symbol)] for i, symbol in enumerate(text)]
#     Pos.sort(key = operator.itemgetter(1))
#     (Pos, Prm) = (list(i) for i in list(zip(*Pos)))

#     # init boolean arr
#     Bh, B2h = ([False for _ in range(n)] for _ in range(2))
#     Bh.append(True)

#     # set degenerate ranks for suffixes with same first chars
#     temp = Prm[0]       
#     Prm[0] = 0
#     Bh[0] = True        
#     for i in range(1, n):
#         if Prm[i] == temp:
#             Prm[i] = Prm[i-1]
#         else:
#             temp = Prm[i]
#             Prm[i] = i
#             Bh[i] = True

#     #init count array
#     Count = [0 for _ in range(n)]
#     # init setup done


#     # while loop
#     h = 1

#     while h<n:
#         # update Prm <- inverse(Pos)
#         for i in range(n):
#             Prm[Pos[i]] = i

#         #reset count array
#         Count = [0 for _ in range(n)]

#         # set Prm[d] <- Prm[head of bucket[d]]
#         head = 0
#         for i in range(1, n):
#             if Bh[i]:
#                 head = Prm[Pos[i]]
#             else:
#                 Prm[Pos[i]] = head   

#         ##  update the suffix that ends this round... ? 
#         # not sure if this is actually? doing? anything?
#         d = n-h
#         e = Prm[d]
#         Prm[d] = e + Count[e]
#         B2h[Prm[d]] = True

#         # For each '=h bucket', update paired keys, set Bh for new heads that arise...
#         left = 0
#         c = 0
#         while c < n:
#             d = Pos[c] - h
#             if d >= 0:
#                 e = Prm[d]
#                 Prm[d] = e + Count[e]
#                 Count[e] += 1
#                 B2h[Prm[d]] = True
#             c += 1
#             # Bh[c] ->if True, c is end of h-bucket; fix B2h
#             if Bh[c]: 
#                 for i in range(left,c):
#                     d = Pos[i]-h
#                     # B2h(Prm[d]) is a new head, 
#                     if d >= 0 and B2h[Prm[d]]: 
#                         j = Prm[d] + 1
#                         # set others in bucket below d to false, end at next head (Bh[j]=Tr)
#                         while True:
#                             if Bh[j] or not B2h[j]:
#                                 break
#                             B2h[j] = False
#                             j += 1
#             # set the left border for the new bucket
#                 left = c

#         # Pos vector <- updated from Prm
#         for i in range(n):
#             Pos[Prm[i]] = i
#         # Bh <- B2h
#         counter=0
#         for i in range(n):
#             if B2h[i]:
#                 Bh[i] = True
#             if Bh[i]:
#                 counter += 1
#         if counter == n:
#             print('\nFinished early, at stage h=',h, 'time = ', time.perf_counter()-start)
#             return Pos
        
#         h = 2*h

#     print('done at h=',h)   
#     return Pos



#####

with open("/Users/jasonmoggridge/Dropbox/Rosalind/Coursera_textbook_track/Course6/data/SuffixArray.txt",'r') as infile:
    text = infile.readline().strip()
#     solutions = list(int(i) for i in infile.readline().strip().split(', '))
    
text_bwt = BurrowsWheeler(text)
# # text_sa = manber_myers(text)
