#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 18:13:03 2019
@author: jasonmoggridge

BA4_H:
    
    SPECTRAL CONVOLUTION
"""



def Spectral_Convolution(spectrum):

    diffs = set()
    multiplicities = { }
    for i in range(len(spectrum)-1):
        for j in range(i+1,len(spectrum)):
            
            diff = abs(spectrum[j]-spectrum[i])
            if diff != 0:
                diffs.add(diff)
                
                if diff in multiplicities:
                    multiplicities[diff] += 1
                else:
                    multiplicities[diff] = 1

    print(multiplicities)
    convolution = []
    print(max(multiplicities.values()))
    while multiplicities:
        
        multi_most = int(max(multiplicities.values()))
        for mass in diffs:
            if mass in multiplicities:
                if multiplicities[mass] == multi_most:
                    for _ in range(multiplicities[mass]):
                        convolution.append(mass)
                    multiplicities.pop(mass)
        
    return convolution

#
""""""
#

#with open("/Users/jasonmoggridge/Desktop/rosalind_ba4h.txt",'r') as f:
#    spectrum = sorted(list(int(i) for i in f.readline().split(' ')))
#    f.close()


spectrum = [0, 86, 160, 234, 308, 320, 382]
convolution = Spectral_Convolution(spectrum)


with open("/Users/jasonmoggridge/Desktop/rosalind_ba4h_output.txt",'w') as out:
    for mass in convolution:
        out.write(str(mass)+' ')
    
    out.close()
