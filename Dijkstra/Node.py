## Dijkstra algorithm realization and variations
# Author: Yang Yu
# Date: 2018/8/9
# Description: Node definition

import math
from math import pi
import numpy as np

## Functions
# Routing nodes function
class Node:
    NodeID = -1
    NodeNum = 0 # Number of all nodes within the map
    Nodepos = np.zeros(2)
    NodeConnections = np.zeros(NodeNum) # Connections from the current node with other nodes (A vec is used to contain the weightings/costs) # All np.inf
    
    # Node internal characteristics
    NodeRoutes = np.zeros([NodeNum,NodeNum]) # Defined as from one outsider node to another. Cost function
    NodePassval = 0
    
    def __init__(self,ID=-1,Pos=[],NodeNum=0,NodePassval=1): # node initialized as passable
        self.NodeID = ID
        self.Nodepos = Pos
        self.NodeNum = NodeNum
        self.NodeConnections = np.ones(NodeNum)*np.inf # Define the internal node structure
        self.NodeRoutes = np.zeros([NodeNum,NodeNum])
        self.NodePassval = NodePassval
        
    def copy(self):
        Node_copy = Node()
        Node_copy.NodeID = self.NodeID
        Node_copy.NodeNum = self.NodeNum
        Node_copy.Nodepos = self.Nodepos.copy()
        Node_copy.NodeConnections = self.NodeConnections.copy()
        Node_copy.NodeRoutes = self.NodeRoutes.copy()
        Node_copy.NodePassval = self.NodePassval
        
        return Node_copy