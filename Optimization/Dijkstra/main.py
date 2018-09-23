## Global Routing algorithm main
# Author: Yang Yu
# Date: 2018/5/2
# Description: NA

import math
from math import pi
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage

import time

from tempfile import TemporaryFile

import matplotlib.pyplot as plt

from Dijkstra import *
import disttable as DT

## Functions
#

## Configurations
MapM = 100
MapN = 100
NodeNum = 100 # Number of nodes
newMapflag = 0

radiusTH = 30

## Form the simulation network (An image based approach)
tmpfile = 'tmp.npz'
if newMapflag == 1: # Create new map
    Map = np.zeros([MapM,MapN])
    NodesX = np.random.randint(0,MapN,NodeNum)
    NodesY = np.random.randint(0,MapM,NodeNum) # Maximum 300 nodes
    Nodespos = np.zeros([2,NodeNum])
    Nodespos[0,:] = NodesX.copy()
    Nodespos[1,:] = NodesY.copy() # Generate the nodes
    
    Nodespos = DT.indx_sort(Nodespos) # Sorting the nodes
    Map = DT.dispMap(Nodespos,MapN,MapM) # Set the nodes to 1
    
    # Generate the edge map (each contains a weighting)
    # This is done here by assigning the distance as the weighting
    # In reality the weighting shall be a lot different
    EdgesW = DT.dist_table_gen(Nodespos)
        
    # Disconnet all edges that are larger than a radius threshold
    indx = np.where(EdgesW > radiusTH)
    EdgesW[indx] = np.inf
    
    # generate the node routes costs (for a given node from one direction to another)
    Nodes = []
    for nodeindx in range(NodeNum):
        NodeID = nodeindx
        Nodepos = Nodespos[:,nodeindx]
        NodeConnections = EdgesW[nodeindx,:]
        
        node = Node(NodeID,Nodepos,NodeNum)
        node.NodeConnections = NodeConnections.copy()
        Nodes.append(node)
        
    # Temporary save
    np.savez(tmpfile,Nodespos=Nodespos,EdgesW=EdgesW,Nodes=Nodes,Map=Map)
else:
    # Temporary load
    npzfile = np.load(tmpfile)
    
    Nodespos = npzfile[npzfile.files[0]]
    EdgesW = npzfile[npzfile.files[1]]
    Nodes = npzfile[npzfile.files[2]]
    Map = npzfile[npzfile.files[3]]
    
plt.figure()
plt.plot(Nodespos[0,:],Nodespos[1,:],'*r')
plt.figure()
plt.imshow(ndimage.interpolation.rotate(Map,90))

## Dijkstra algorithm simu
StrNodeIndx = 0
EndNodeIndx = len(Nodes)-1 # The algorithm would try to find the minimum distance route to the endpt

# Nodes_dist,Nodes_prev = Dijkstra_Path_Dist(Nodes,StrNodeIndx)
path,path_length,pathpos = Dijkstra(Nodes,StrNodeIndx,EndNodeIndx)
    
plt.figure()
plt.plot(Nodespos[0,:],Nodespos[1,:],'*r')
plt.plot(pathpos[0,:],pathpos[1,:],'--')

##
plt.show()
