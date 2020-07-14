#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 22:56:56 2019
@author: jasonmoggridge

http://rosalind.info/problems/ba3g/

BA3_G:
    Eularian path in digraph
    
    
"""


import sys
sys.setrecursionlimit(10**6) 

"""Find the source node to initiate dfs at:"""
def get_source(Adj):
    deg = {}
    
    for v in Adj:
        
        if v not in deg.keys():
            deg[v] = -len(Adj[v])
        else:
            deg[v] -= len(Adj[v])
    
        for u in Adj[v]:
            if u not in deg.keys():
                deg[u] = 1
            else:
                deg[u] += 1
    for v in deg:
        if deg[v] == -1:
            source = v 
    return source


"""Eularian path/circuit dfs"""
def dfs(v, graph, path, stack):
    
    print('dfs: v:',v,'\nPath:',path,'\nStack:',stack)

    if v not in graph.keys():
        print('\thit the sink')
        path.append(v)
        if len(stack) == 0:
            print('\t\tNo stack left\n\t\tPath is:', path)
            return
        else:
            print('\t\tStill stack >>>', stack[-1])
            dfs(stack.pop(), graph, path, stack)
        
    elif len(graph[v]) == 0:
        print('\tNo neighbors')
        path.append(v)
        if len(stack) == 0:
            print('\t\tNo stack left\n\t\tPath is:', path)
            return
        else:
            print('\t\tStill stack >>>', stack[-1])
            dfs(stack.pop(), graph, path, stack)
    else:
        print('\t Going deeper, appending stack, going to :', graph[v][0])
        stack.append(v)
        dfs(graph[v].pop(0), graph, path, stack)
    ###
    



f = open('//Users/jasonmoggridge/Desktop/rosalind_ba3g.txt', 'r')
lines = list(str(l.strip('\n')) for l in f.readlines())
Edges = [[i.strip(' ') for i in line.split('->')] for line in lines]

Adj = {}
for edge in Edges:
    Adj[int(edge[0])]=[int(i) for i in edge[1].split(',')]
del((Edges, lines, f))


source = get_source(Adj)
stack = []
path = []
dfs(source, Adj, path, stack)

path.reverse()
string = '->'.join(str(i) for i in path)
#print('\n\n\n',string)

with open('//Users/jasonmoggridge/Desktop/rosalind_ba3g_output.txt', 'w') as wr:
    wr.write(string)
    wr.close()
