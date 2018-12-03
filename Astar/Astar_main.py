## Global Routing algorithm main (Astar Heuristic)
# Author: Yang Yu
# Date: 2018/5/2
# Description: A star algorithm

import math
from math import pi
import numpy as np
import scipy as sp
import scipy.ndimage as SPDI
import time

from tempfile import TemporaryFile
import matplotlib.pyplot as plt
import meshop as MP

## Functions
class Node: #
    
    NodeID = -1 # Node ID
    x = -1 # x position
    y = -1 # y position
    
    # Connections
    ConnectionsID = [] # PrevConnections nods (NodeID)
    Connections_cost= [] # Connections cost, travel from a node to another
    Connections_dir = [] # To record the directions from this node to another
    
    # Metrics
    prevNode = -1
    hmetrics = np.inf
    gmetrics = np.inf
    fmetrics = np.inf
    
    def __init__(self,NodeID,x=-1,y=-1,ConnectionsID=[],Connections_cost=[],Connections_dir=[]):
        self.NodeID = NodeID
        self.x = x
        self.y = y
        self.ConnectionsID = ConnectionsID.copy()
        self.Connections_cost = Connections_cost.copy()
        self.Connections_dir = Connections_dir.copy()
    
    def copy(self):
        NodeCopy = Node(self.NodeID)
        NodeCopy.x = self.x
        NodeCopy.y = self.y
        NodeCopy.ConnectionsID = self.ConnectionsID.copy()
        NodeCopy.Connections_cost = self.Connections_cost.copy()
        NodeCopy.Connections_dir = self.Connections_dir.copy()
        
        return NodeCopy
    
    def update_prevNode(self,prevNode):
        self.prevNode = prevNode
    
    def metrics_update(self,gmetrics,hmetrics,fmetrics):
        self.gmetrics = gmetrics
        self.hmetrics = hmetrics
        self.fmetrics = fmetrics
    
def GMap2Nodes(GMap,NodeIDmap,ConnectionTemplate):
    
    M,N = GMap.shape
    
    GNodes = []
    NodeID = -1
    for xx in range(M):
        for yy in range(N):
            
            # print(xx,yy)
            
            if GMap[xx,yy] == 0:
                continue
            elif GMap[xx,yy] == 1:
                NodeID = NodeID+1
                newNode = Node(NodeID,xx,yy)
                
                # Determine the connections
                for ii in range(3):
                    for jj in range(3):
                        
                        if ii == 1 and jj == 1:
                            continue
                        
                        xx_n = xx-1+ii
                        yy_n = yy-1+jj # neighbour candidate
                        
                        if xx_n < 0 or xx_n > M-1 or yy_n < 0 or yy_n > N-1:
                            continue
                        elif GMap[xx_n,yy_n] == 0:
                            continue
                        else:
                            newNode.ConnectionsID.append(np.uint(NodeIDmap[xx_n,yy_n]))
                            newNode.Connections_cost.append(ConnectionTemplate[ii,jj])
                            newNode.Connections_dir.append([ii-1,jj-1])
                
            GNodes.append(newNode.copy())
            
    return GNodes
    
def angle_alter_cost_cal(relative_dir1,relative_dir2,ang_alt_cost_pi = 0):
    
    ang1 = np.arctan2(-relative_dir1[1],-relative_dir1[0])
    ang2 = np.arctan2(relative_dir2[1],relative_dir2[0]) # Only with regard to direction
    ang_alt = ang2-ang1
    
    # ang_alt to cost relation
    ang_alt_cost = ang_alt_cost_pi/(pi)*ang_alt
    
    return ang_alt_cost

def imrotate(I,ang):
    
    Ir = SPDI.interpolation.rotate(I,ang,reshape=0,order=1)
    
    return Ir

def vec_add_element(vec,a):
    
    vec_l = vec.tolist()
    vec_l.append(a)
    vec_e = np.array(vec_l)
    
    return vec_e
    
def list_add(l1,l2):
    
    l1v = np.array(l1)
    l2v = np.array(l2)
    l12v = l1v+l2v
    l12 = l12v.tolist()
    
    return l12
    
# A-star search algorithm
# This is a modified nodes version
def Astar(GNodes,StrNodeIndx,EndNodeIndx):
    
    NodesNum = len(GNodes)
    
    # Closed set: nodes evaluated set
    closedsetID = []
    
    # Openset: nodes to be evaluated
    opensetID = []
    gmetrics_pool = [] # Travel distance: travel metrics includes the turning and the traveling cost
    hmetrics_pool = [] # Heuristic estimation of the openset
    fmetrics_pool = [] # Full metrics: combination of hmetrics and gmetrics
    
    GNode_str = GNodes[StrNodeIndx]
    GNode_end = GNodes[EndNodeIndx]
    
    # Include the first element in the openset
    opensetID.append(GNode_str.NodeID)
    gmetrics_pool.append(0)
    hmetrics_pool.append(hmetrics_estimate(GNodes,GNode_str.NodeID,EndNodeIndx))
    fmetrics_pool.append(gmetrics_pool[-1]+hmetrics_pool[-1])
    
    # Update the GNode's metrics
    GNodes[StrNodeIndx].metrics_update(gmetrics_pool[-1],hmetrics_pool[-1],fmetrics_pool[-1])
    
    # Routing
    success_indx = 0
    while len(opensetID) > 0: # openset not empty
        
        print(opensetID)
        
        # Search the openset for the best candidate
        fmetrics_indx = np.argmin(np.array(fmetrics_pool)) # Find the best element in the openset
        current_candidate = GNodes[opensetID[fmetrics_indx]] # Current candidate
        current_gmetrics = gmetrics_pool[fmetrics_indx]
        current_hmetrics = hmetrics_pool[fmetrics_indx]
        current_fmetrics = fmetrics_pool[fmetrics_indx]
        
        print(current_candidate.x,current_candidate.y)
        print(current_candidate.NodeID)
        
        # Check if the goal has been reached
        if current_candidate.NodeID == GNode_end.NodeID:
            success_indx = 1
            # A route could be established by tracing from current candidate's prevNode
            break # Goal reached
        else:
            # Update the closed set
            closedsetID.append(current_candidate.NodeID)
            
            # Remove the ID from opensetID (and corresponding metrics)
            opensetID.pop(fmetrics_indx)
            gmetrics_pool.pop(fmetrics_indx)
            hmetrics_pool.pop(fmetrics_indx)
            fmetrics_pool.pop(fmetrics_indx)
        
        print(closedsetID)
        
        # Update the openset
        openset_candidates_ID = current_candidate.ConnectionsID.copy()
        for ii in range(len(openset_candidates_ID)):
            
            openset_candidate_ID = openset_candidates_ID[ii]
            openset_candidate = GNodes[openset_candidate_ID] # Pointer for convenience
            
            # Check if the openset candidate is in the closed set
            InCloseSetFlag,_ = CheckIfInSet(closedsetID,openset_candidate_ID)
            
            if InCloseSetFlag == 0: # Empty: not in closed set
                NodeID_current = current_candidate.NodeID
                NodeID_prev = current_candidate.prevNode
                NodeID_next = openset_candidate.NodeID
                # To calculate the value from the prevNode through the current node to the openset candidate
                
                # Calculate the gmetrics and hmetrics (and fmetrics)
                cand_gmetrics = gmetrics_estimate(GNodes,NodeID_prev,NodeID_current,NodeID_next)
                cand_hmetrics = hmetrics_estimate(GNodes,NodeID_next,EndNodeIndx)
                cand_fmetrics = cand_gmetrics + cand_hmetrics
                
                # Check if the openset candidate is in the closed set
                InOpenSetFlag,InOpenSetIndex = CheckIfInSet(opensetID,openset_candidate_ID)
                
                if InOpenSetFlag == 0: # not in open set
                    # Add the openset candidate into the openset
                    # Update the metrics
                    openset_candidate.update_prevNode(current_candidate.NodeID)
                    openset_candidate.metrics_update(cand_gmetrics,cand_hmetrics,cand_fmetrics)
                    opensetID.append(openset_candidate.NodeID)
                    gmetrics_pool.append(cand_gmetrics)
                    hmetrics_pool.append(cand_hmetrics)
                    fmetrics_pool.append(cand_fmetrics)
                else: # in the open set
                    if cand_gmetrics < openset_candidate.gmetrics:
                        openset_candidate.update_prevNode(current_candidate.NodeID)
                        openset_candidate.metrics_update(cand_gmetrics,cand_hmetrics,cand_fmetrics)
                        gmetrics_pool[InOpenSetIndex[0][0]] = cand_gmetrics # update the metrics
                        hmetrics_pool[InOpenSetIndex[0][0]] = cand_hmetrics
                        fmetrics_pool[InOpenSetIndex[0][0]] = cand_fmetrics
                
    return success_indx,current_candidate
    
def Astar_routing(GNodes,current_candidate):
    
    route = []
    route.append(current_candidate.NodeID) # Adding the final Node ID
    
    candidate = current_candidate
    while 1:
        prevNode = candidate.prevNode
        route.append(prevNode)
        
        if prevNode == -1:
            break
        else:
            candidate = GNodes[prevNode]
    
    return route
    
def gmetrics_estimate(GNodes,NodeID_prev,NodeID_current,NodeID_next,init_dir = np.array([1,1])):
    
    Nodecurrent = GNodes[NodeID_current]
    Nodenext = GNodes[NodeID_next]
    if NodeID_prev == -1:
        Nodeprev = GNodes[NodeID_current]
    else:
        Nodeprev = GNodes[NodeID_prev]

    xcurrent = Nodecurrent.x
    ycurrent = Nodecurrent.y
    xnext = Nodenext.x
    ynext = Nodenext.y
    
    # Calculate the distance between previous and current nodes
    dist_current_next = np.sqrt((xnext-xcurrent)**2+(ynext-ycurrent)**2)
    gmetrics_tentative = dist_current_next
    
    # Rotation
    CurrentConnectionsID = Nodecurrent.ConnectionsID
    # Find the direction from current to next
    indx = np.where(np.array(CurrentConnectionsID)==Nodenext.NodeID)
    CurrentNextConnectiondir = np.array(Nodecurrent.Connections_dir[indx[0][0]])
    # Find the direction from current to prev
    if NodeID_prev == -1:
        PrevCurrentConnectiondir = init_dir
    else:
        indx = np.where(np.array(CurrentConnectionsID)==Nodeprev.NodeID)
        PrevCurrentConnectiondir = np.array(Nodecurrent.Connections_dir[indx[0][0]])
    
    ang_rot_cost = angle_alter_cost_cal(PrevCurrentConnectiondir,CurrentNextConnectiondir)
    
    gmetrics = GNodes[NodeID_current].gmetrics+gmetrics_tentative+ang_rot_cost
    
    return gmetrics
    
def hmetrics_estimate(GNodes,CurrentNodeID,EndNodeIndx,effect=1):
    
    if effect == 0:
        hmetrics = 0
    else:
        # Find the distance between the NodeID and the endNode
        CurrentNode = GNodes[CurrentNodeID]
        EndNode = GNodes[EndNodeIndx]
        
        xcurrent = CurrentNode.x
        ycurrent = CurrentNode.y
        xend = EndNode.x
        yend = EndNode.y
        
        dist_current_end = np.sqrt((xend-xcurrent)**2+(yend-ycurrent)**2)
        hmetrics = dist_current_end
        
    return hmetrics

def CheckIfInSet(closedsetID,openset_candidate_ID):
    
    closedsetID_A = np.array(closedsetID)
    indx_l = np.abs(closedsetID_A-openset_candidate_ID)
    index = np.where(indx_l == 0)
    
    flag = 1
    if len(index[0]) == 0:
        flag = 0
    
    return flag,index

def plot_rect(x_range,y_range):
    
    line1 = np.array([[x_range[0],x_range[0]],[y_range[0],y_range[-1]]])
    line2 = np.array([[x_range[0],x_range[-1]],[y_range[0],y_range[0]]])
    line3 = np.array([[x_range[-1],x_range[-1]],[y_range[0],y_range[-1]]])
    line4 = np.array([[x_range[0],x_range[-1]],[y_range[-1],y_range[-1]]])
    
    plt.plot(line1[0],line1[1])
    plt.plot(line2[0],line2[1])
    plt.plot(line3[0],line3[1])
    plt.plot(line4[0],line4[1])

## Configurations
MapM = 50
MapN = 50

# In this case we form a map that contains all pixel
newMapflag = 1

# Get the global map
GMap = np.ones([MapM,MapN]) # M y-axis, N x-axis (x-axis from upper to down, y from left to right)
# Rotate 90 degree to align with x,y plots
indx = np.where(GMap == 0)
GMap[indx] = -1

# Obstacle definition
x_range = np.arange(15,35)
y_range = np.arange(15,25)
for x in x_range:
    for y in y_range:
        GMap[x,y] = 0

NodeNum = MapM*MapN-len(indx[0]) # Calculate the possible nodes' number

# Connection template generation
ConnectionTemplate = np.zeros([3,3])
mg = MP.meshgridi(ConnectionTemplate,C=[1,1])
ConnectionTemplate = mg.meshr.copy()*10 # Distance is multiplied by 10

NodeIDmap = np.zeros(GMap.shape)
NodeID = 0
for xx in range(MapM):
    for yy in range(MapN):
        
        if GMap[xx,yy] == 0:
            continue
        elif GMap[xx,yy] == 1:
            NodeID = NodeID+1
            NodeIDmap[xx,yy] = NodeID
NodeIDmap = NodeIDmap-1

# plt.figure()
# plt.imshow(imrotate(NodeIDmap,ang))

## Form the simulation network (An image based approach)
tmpfile = 'Astar_tmp.npz'
if newMapflag == 1: # Create new map
    
    GNodes = GMap2Nodes(GMap,NodeIDmap,ConnectionTemplate)
    
    # Temporary save
    np.savez(tmpfile,GNodes=GNodes)
else:
    # Temporary load
    npzfile = np.load(tmpfile)
    
    GNodes = npzfile[npzfile.files[0]]

# Test of angle rotation cost
# relative_dir1 = np.array([0,1])
# relative_dir2 = np.array([0,1])
# angle_alter_cost = angle_alter_cost_cal(relative_dir1,relative_dir2)

## Dijkstra algorithm simu
StrNodeIndx = 0
EndNodeIndx = len(GNodes)-1 # The algorithm would try to find the minimum distance route to the endpt
# From one corner to another corner

success_indx,current_candidate = Astar(GNodes,StrNodeIndx,EndNodeIndx)
route = Astar_routing(GNodes,current_candidate)

route_x = []
route_y = []
indx = len(route)-1
for ii in range(len(route)-1):
    route_x.append(GNodes[route[indx-ii-1]].x)
    route_y.append(GNodes[route[indx-ii-1]].y)

plt.figure()
plot_rect(x_range,y_range)
plt.plot(route_x,route_y,'r*')

##
plt.show()
