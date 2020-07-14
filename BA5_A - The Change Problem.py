#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 01:30:41 2020

@author: jasonmoggridge

Change Problem: Find the minimum number of coins needed to make change.

Input: An integer money and an array Coins of d positive integers.
Output: The minimum number of coins with denominations Coins that changes money.



    GreedyChange(money)
        change ← empty collection of coins)
        while money > 0
            coin ← largest denomination that is less than or equal to money
            add a coin with denomination coin to the collection of coins change
            money ← money − coin
        return change
        
   
    RecursiveChange(money, Coins)
        if money = 0
            return 0
        MinNumCoins ← ∞
        for i ← 0 to |Coins| - 1
            if money ≥ coini
                NumCoins ← RecursiveChange(money − coini, Coins)
                if NumCoins + 1 < MinNumCoins
                    MinNumCoins ← NumCoins + 1
         return MinNumCoins
        
        

MakeChange - dynamic prog

    This algo is memoized but not efficient bc it tries to search through the
    entire tree in a top down fashion. 
    
    It is far more efficient to build up the array from 0 to money
    Rather than money -> zero by trying the smallest denomination in each step
    
    
    19998
21,17,10,5,3,1

 
"""
#import sys
#sys.setrecursionlimit(10**7)
#
#money = 76
#coins =[5, 4, 1]
##change = [int(0) for i in range(len(coins))]
#memo = {0:0}
#
#
#def MinNumCoins(money):
#    
#    if money in memo:
#        return memo[money]
#
#    else:
#        paths = []
#        for coin in coins:
#            if money - coin >= 0:
#                paths.append(MinNumCoins(money-coin)+1)
#        memo[money] = min(paths)
#        return memo[money]  
#    
#coiners = MinNumCoins(money)


"""
   DPChange(money, Coins):
     MinNumCoins(0) ← 0
     for m ← 1 to money
        MinNumCoins(m) ← ∞
            for i ← 0 to |Coins| - 1
                if m ≥ coini
                    if MinNumCoins(m - coini) + 1 < MinNumCoins(m)
                        MinNumCoins(m) ← MinNumCoins(m - coini) + 1
    output MinNumCoins(money)
    
"""

def DPChange(money, coins):
    
    minCoins = {0:0}
    for m in range(1,money+1):
        minCoins[m] = float('inf')
        for coin in coins:
            if m >= coin:
                if minCoins[m-coin]+1 < minCoins[m]:
                    minCoins[m] = minCoins[m-coin]+1
    print(minCoins)
    return minCoins[money]
##

import sys
sys.setrecursionlimit(10**7)
money = 19998
coins =[21,17,10,5,3,1]
dp = DPChange(money, coins)




"""
Recall that our original goal was to make change, not just compute
 MinNumCoins(money). 
 Modify DPChange so that it not only computes
 the minimum number of coins but also returns these coins.

Having problems with 2d list where trying to update list caused all
lists to change, including the zero value
    
-> any time you do:
            listY = listX
            listY[i] = p

    -> listX[i] will also change bc listY is just a pointer to list X,
    to simply copy the values in the list, use:
            listY = copy.deepcopy(listX)
        and then listX won't change when listY is edited

"""

# this took forever to figure out. never forget how copying works in lists ok
# have to use deepcopy(list) to get its values without lists all changing together.

def DPChangeArray(money, coins):

    import copy
    minCoins = {0:[0 for coin in coins]}
    minCoins.update(zip(range(1,money+1),\
                        ([float('inf') for coin in coins] for i in range(1,money+1))))
    for m in range(1,money+1):        
        for coin in coins:
            if m >= coin:
                if sum(minCoins[m - coin]) + 1 < sum(minCoins[m]):
                    minCoins[m] = copy.deepcopy(minCoins[m - coin])
                    minCoins[m][coins.index(coin)] += 1
    return minCoins[m]
###
money = 69
coins =[21,17,10,5,3,1]
minCoins = DPChangeArray(money, coins)





"""
STOP and Think: Consider the following modifications of DPChange.

If money = 10^9, DPChange requires a huge array of size 10^9. 
Modify the DPChange algorithm so that the array size required
 does not exceed the value of the largest coin denomination.
 """
 
 # will solve later, getting really sick of the make change problem. 

#



