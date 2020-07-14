#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 17:05:28 2019

@author: jasonmoggridge

BA4-exercise:
    Exercise Break: How many subpeptides does a linear peptide
    of given length n have? 
    (Include the empty peptide and the entire peptide.)

Input: An integer n.
Output: The number of subpeptides of a linear peptide of length n.



Number of subpeptides for a length n linear peptide:
    is the sum of range(1,2,....n) => n*(n+1)/2 , add first element empty string

   Number = n*(n+1)  +  1
            ------
               2
   
   
"""

# length n
n = 100

linear_subpeptides = n*(n+1)//2 + 1
