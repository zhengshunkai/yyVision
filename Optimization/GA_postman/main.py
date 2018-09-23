## GA postman problem
# Author: Yang Yu
# Date: 2018/9/23
# Description: a GA postman problem solution based on GA algorithm

import math
from math import pi
import numpy as np
import scipy as sp
import scipy.ndimage as ndimage
import time
from tempfile import TemporaryFile
import matplotlib.pyplot as plt
import random
from copy import deepcopy
from rearrange import rearrange

import disttable as DT
from Node import Node

## Functions  
def GA_init(Nodes):
    # The solution is essentially a combination of all nodes
    NodeNum = len(Nodes)

    init_combinations = []
    comb0 = np.arange(NodeNum-2) # Always from node 0 to node N
    init_combinations.append(comb0+1) # From node 1 to node N
    
    # Create M random order
    # Get 5 random candidates
    for ii in range(5):
        comb1 = comb0.copy()
        comb1 = shuffle_array(comb1)
        init_combinations.append(comb1.copy()+1)
        
    return init_combinations

# Genre operations here are naive
def GA_selection(Nodes,genre_population,selNum):
    
    genreNum = len(genre_population)
    metrics_pool = np.ones(genreNum)*np.inf
    for ii in range(genreNum):
        metrics_pool[ii] = GA_metrics(Nodes,genre_population[ii])
    
    print(metrics_pool)
    
    # Select the top selNum candidates
    indx = np.argsort(metrics_pool)
    if genreNum <= selNum:
        new_genre_population = rearrange(genre_population,indx)
    else:
        new_genre_population = []
        for ii in range(selNum):
            new_genre_population.append(genre_population[indx[ii]].copy())
    
    return new_genre_population
    
def GA_metrics(Nodes,genre):
    
    # Calculate the overall distance of the seq.
    NodeNum = len(Nodes)
    dist_overall = 0
    Node0 = Nodes[0].copy()
    for ii in range(NodeNum-2):
        Node1 = Nodes[genre[ii]].copy()
        disttmp = np.sqrt(np.sum((Node0.Nodepos-Node1.Nodepos)**2))
        dist_overall = dist_overall+disttmp
        Node0 = Node1.copy()
    
    Node1 = Nodes[-1].copy()
    disttmp = np.sqrt(np.sum((Node0.Nodepos-Node1.Nodepos)**2))
    dist_overall = dist_overall+disttmp
    
    return dist_overall
    
def GA_mutate(Nodes,genre):    
    
    # Switch 2 position in the genre
    NodeNum = len(Nodes)
    pos1 = random.randint(0,NodeNum-3)
    while 1:
        pos2 = random.randint(0,NodeNum-3)
        if pos1 != pos2:
            break
        print(pos1,pos2)
    
    genre_mut = genre.copy()
    genre_mut[pos1],genre_mut[pos2] = genre_mut[pos2],genre_mut[pos1]
    
    return genre_mut

def GA_mate(Nodes,genre0,genre1):
    
    genre_mat = genre0.copy()
    NodeNum = len(Nodes)
    
    genre_sum = genre0 + genre1
    genre_mat = np.argsort(genre_sum)+1
    
    return genre_mat

def GA_routing(Nodes,itermax=100,selNum=4,mateNum=2,mutNum=3):
    
    # Initialization
    prev_population = GA_init(Nodes)
    for iter in range(itermax):
        # Selection
        population = GA_selection(Nodes,prev_population,selNum)
        
        # Mating
        for ii in range(mateNum):
            female = population[ii]
            malepos = np.random.randint(ii+1,selNum)
            male = population[malepos]
            
            candidate_mate = GA_mate(Nodes,female,male)
            population.append(candidate_mate.copy())
        
        # Mutation
        for ii in range(mutNum):
            candidate_mut = GA_mutate(Nodes,population[ii])
            population.append(candidate_mate)
        
        prev_population = deepcopy(population)
        
        # Early termination
        # TBD
        
    population = GA_selection(Nodes,prev_population,selNum)
    GA_genre = population[0]
    
    return GA_genre,population
    
def shuffle_array(x):
    
    xl = x.tolist()
    random.shuffle(xl)
    x = np.array(xl)
    
    return x

## Configurations
MapM = 100
MapN = 100
NodeNum = 10 # Number of nodes
newMapflag = 1 # New map stored in tmp.npz

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
    # For postman the last point overlap with the first
    Nodespos[:,0] = Nodespos[:,-1]
    
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
plt.plot(Nodespos[0,:],Nodespos[1,:],'b--')
plt.figure()
plt.imshow(ndimage.interpolation.rotate(Map,90))

# GA routines
final_combination,final_population = GA_routing(Nodes,selNum=4)
Nodespos_final = np.zeros([2,NodeNum])
for ii in range(NodeNum):
    if ii == 0 or ii == NodeNum-1:
        Nodespos_final[:,ii] = Nodes[ii].Nodepos 
    else:
        Nodespos_final[:,ii] = Nodes[final_combination[ii-1]].Nodepos 

plt.figure()
plt.plot(Nodespos[0,:],Nodespos[1,:],'*r')
plt.plot(Nodespos_final[0,:],Nodespos_final[1,:],'b--')

#
plt.show()