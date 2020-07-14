#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 21:36:20 2020
@author: jasonmoggridge


BA-7a:
Compute Distances Between Leaves 

In this chapter, we define the length of a path in a tree as the sum
 of the lengths of its edges (rather than the number of edges on the path).
 As a result, the evolutionary distance between two present-day species
 corresponding to leaves i and j in a tree T is equal to the length
 of the unique path connecting i and j, denoted di, j(T).

Distance Between Leaves Problem
Compute the distances between leaves in a weighted tree.

Given: An integer n followed by the adjacency list of a weighted tree
 with n leaves.

Return: A space-separated n x n (di, j), where di, j is the length of
 the path between leaves i and j.

Sample Dataset

4
0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6


Sample Output

0   13  21  22
13  0   12  13
21  12  0   13
22  13  13  0


IMPLEMENT DIJKSTRA'S ALGORITHM
    # other option is to use Floyd-Warshal algorithm if graph
    # is highly connected/ ie. not sparse like a tree


"""


def Dijkstra(adjacency, source): 
    """
    Dijkstra takes adjacency list, source node
    returns an array of distances of shortest paths source -> i
    """

    def Decrease_key(heap, v):   
    
        """
        Maintains the order of nodes in min_heap after an edge is visited and
        the distance(source->node) is updated. Node is moved to correct position in 
        min heap queue.      
    
        """
        # Find the node and remove it from heap
        # items in heap are (node, Distance[node])
        for node in heap:
            if node[0] == v:
                heap.remove(node)
                break

        for node in heap:
            if distance[node[0]] > distance[v]:
                ins = heap.index(node)
                heap.insert(ins, (v, distance[v]))
                break


    # initial distance array ->[ 0 for source, inf for rest of nodes]
    distance = [float('inf') for _ in range(len(adjacency))]
    distance[source] = 0
    
    # intial heap -> source as first node to visit
    # Implement this ..(if graph is huge esp.,) can start with heap = (source,0)
    # and just add (node, dist) as they are first seen in BFS
    
    heap = [(node, distance[node]) for node in range(len(adjacency))]            
    heap.remove((source,0))
    heap = [(source, 0)] + heap

    
    # BFS until all edges are relaxed by visiting all nodes in heap (while heap:)
        # Pop next node (min current dist) from top of heap
        # Visit all outgoing edges and update distance[u] if path source->v->u is shorter
        # if source -> u is relaxed, decr_key(u)
    
    while heap:
        v = heap.pop(0)                       
        for u in adjacency[v[0]]:
            if distance[u[0]] > distance[v[0]] + u[1]:
                distance[u[0]] = distance[v[0]] + u[1]
                Decrease_key(heap, u[0])
    return distance 
##

    
# Main 
    

# Parse weighted edges to adjacency list:
    # First line of input is number of leaves
    # nodes labelled from (0...n-1) are leaves (deg=1)
    # nodes (n...) are internal nodes
    # edges input format 'v->u:x' v- vertice, u- neighbor, x- weight

edges= []
with open("/Users/jasonmoggridge/Desktop/rosalind_ba7a.txt",'r') as infile:
    n = int(infile.readline().strip())
    for line in infile.readlines():
        edge = line.split('->')
        dest = edge[1].split(':')
        edges.append(tuple(int(x) for x in (edge[0], dest[0], dest[1] ) ))
    infile.close()

adjacency = {}
for edge in edges:
    if edge[0] not in adjacency.keys():
        adjacency[edge[0]] = [(edge[1:])]
    else:
        adjacency[edge[0]].append((edge[1:]))
del(dest,edge,edges,line)    


# Dijkstra:

# Initialize distance matrix D    
D = []

# call Dijkstra for each node in graph to get shortest paths
# add array as row to matrix
for i in range(len(adjacency)):
    d = Dijkstra(adjacency, i)
    D.append(d)
 
# Output the inner matrix of just leaves -> D[:n][:n]
with open("/Users/jasonmoggridge/Desktop/dataset_output_10328_12.txt",'w') as outfile:
    for i in range(n):
        line = ' '.join(str(j) for j in D[i][:n])
        outfile.write(line + '\n')
        print(line)
outfile.close()
del(d,i,n,line)