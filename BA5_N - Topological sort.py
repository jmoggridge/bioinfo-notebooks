#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 23:24:02 2020

@author: jasonmoggridge
"""



def adj_list(file):

    edges = []
    for edge in file:
        edges.append(edge.strip().split('->'))    
    adjacency = {}
    for edge in edges:
        v = int(edge[0])
        U = [int(i) for i in edge[1].split(',')]
        if v not in adjacency:
            adjacency[v] = U
        else:
            adjacency[v] += U
        for u in U:
            if u not in adjacency:
                adjacency[u] = []
    return adjacency
#

def Topological_Sort(graph): 
# dfs builds topsort upon ret from expl'd v.
    
    def dfs(v, topo, visited, graph): 
        visited[v] = True
        for u in graph[v]:
            if not visited[u]:
                dfs(u, topo, visited, graph)            
        topo.append() 

    visited ={}
    for v in graph:
        visited[v] = False

    topo = []
    for v in graph:
        if not visited[v]:
            dfs(v, topo, visited, graph)
            
    return topo[::-1] 
    # list needs to be inverted #


with open("/Users/jasonmoggridge/Desktop/rosalind_ba5n.txt",'r') as f:
    graph = adj_list(f)
    
topo = Topological_Sort(graph)
print(topo)