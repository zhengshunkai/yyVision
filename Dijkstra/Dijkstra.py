## Dijkstra algorithm realization and variations
# Author: Yang Yu
# Date: 2018/8/9
# Description: Dijkstra algorithm and its variations (e.g. A*)

import math
from math import pi
import numpy as np
import scipy as sp

import time
import matplotlib.pyplot as plt

import disttable as DT
from Node import Node
from comb_index import *

## Functions
# Baseline Dijkstra algorithm
# This is a modified nodes version
def Dijkstra_Path_Dist(Nodes,StrNodeIndx):
    # Calculate the Dijsktra distance between the StrNode and every nodes
    # Calculate edges weightings only
    
    # Initialization
    NodeQueue = np.ones(len(Nodes)) # mark if the node has been used or not
    Nodes_dist = np.ones(len(Nodes))*np.inf
    Nodes_prev = np.ones(len(Nodes))*np.nan
    
    Nodes_dist[StrNodeIndx] = 0 # starting point
    
    total_unsearched_nodes = len(np.where(NodeQueue == 1)[0])
    while total_unsearched_nodes > 0:
        
        indx1 = np.where(NodeQueue == 1)[0] # Excluding the previous points
        indx2 = np.nanargmin(Nodes_dist[indx1]) # Find the minimum distance
        nodeindx = indx1[indx2]
        
        # Get the neighbours of node_sel
        node = Nodes[nodeindx]
        NodeConnections = node.NodeConnections.copy()
        
        indx = np.where((NodeConnections > 0)*(NodeQueue == 1))[0] # Find the neighbours, excluding previous points
        
        for ii in range(len(indx)):
        
            node_NB_indx = indx[ii] # Find the neighbours
            alt = Nodes_dist[nodeindx]+NodeConnections[node_NB_indx] # Update the distance
            
            if alt < Nodes_dist[node_NB_indx]:
                Nodes_dist[node_NB_indx] = alt
                Nodes_prev[node_NB_indx] = nodeindx
        
        NodeQueue[nodeindx] = -1 # Switch that to -1
        total_unsearched_nodes = len(np.where(NodeQueue == 1)[0])
    
    return Nodes_dist,Nodes_prev
    
def Dijkstra(Nodes,StrNodeIndx,EndNodeIndx): # Not Fib heap version
      
    # Find all shortest path from start node to all nodes
    Nodes_dist,Nodes_prev = Dijkstra_Path_Dist(Nodes,StrNodeIndx)
    
    # Track from the EndNode to the StrNode
    path = []
    Nodeindx = EndNodeIndx
    
    if np.isnan(Nodes_prev[Nodeindx]):
        pass # No path
    else:
        while 1:
            path.append(Nodeindx)
            if Nodeindx == StrNodeIndx:
                break
            else:
                Nodeindx = np.int(Nodes_prev[Nodeindx])
    
    path = np.array(path)
    path_length = Nodes_dist[EndNodeIndx]
    
    # Construct the path
    pathpos = []
    for ii in range(len(path)):
        path_node_tmp = Nodes[path[ii]].Nodepos.copy()
        pathpos.append(path_node_tmp)
        
    pathpos = np.array(pathpos).T
    
    return path,path_length,pathpos