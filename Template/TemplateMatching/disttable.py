# Generation of (global/local) disttable 
# Author: Yang Yu
# Date: 2018/4/10
# Description: The global disttable records the relation between each control points
# Each control points is represented by a line/col combination
# Each point in the map contains the value of the distance between two control points

import math
from math import pi
import numpy as np
import scipy as sp

# Configuration (global)
class globalMap_config:
    global_length = [] # meter
    global_width = [] # meter
    simu_prec = [] # meter (1 cm)
    disp_prec = [] # meter (10 cm)
    reflectorNum = []
    
    def __init__(self):
        self.global_length = 50 # meter
        self.global_width = 50 # meter
        self.simu_prec = 0.001 # meter (0.1 cm)
        self.disp_prec = 0.5 # meter (50 cm)
        self.reflectorNum = 50
        
# Distance metrics
def distmetrics(PtA,PtB):
# Metrics shall be rotation invariant
    m = np.sqrt(np.nansum(np.square(PtA-PtB)))
    
    return m
    
def angmetrics(PtA,PtB):
# metrics of angle
    if PtA[0]-PtB[0] == 0:
        ang = np.nan # Not valid
    else:
        ang = np.array(np.arctan2((PtB[0]-PtA[0]),(PtB[1]-PtA[1])),dtype='float')
        
    return ang

def angadjrange(ang):
    
    inflag = 0
    while inflag == 0:
        if ang > pi:
            ang = ang - 2*pi
        elif ang <= -pi:
            ang = ang + 2*pi
        else:
            inflag = 1
    
    return ang
    
# Sort the control points
def indx_sort(pos_unsorted,GMconfig=[]):
    
    # global_length = GMconfig.global_length
    # Sorting rule(s)
    # pos_unsorted_indx = pos_unsorted[0,:]*global_length+pos_unsorted[1,:] # line by line search
    # pos_unsorted_indx = np.sqrt(np.square(pos_unsorted[0,:])+np.square(pos_unsorted[1,:])) # Radius min to max search
    pos_unsorted_indx = np.abs(pos_unsorted[0,:]*pos_unsorted[1,:])
    # pos_unsorted_indx = np.abs(pos_unsorted[0,:])+np.abs(pos_unsorted[1,:]) # coordinates sum
   
    pos_ord = np.argsort(pos_unsorted_indx)
    pos_sorted = np.array(pos_unsorted[:,pos_ord],dtype='float')
    
    return pos_sorted
# Ideally, the distance between two points in the disttable shall represent the distance from the map (say points 1 and 100 shall be farther away from points 1 and 10)

# Generate position information from a map
def map2pos(Map,GMconfig):
    
    disp_prec = GMconfig.disp_prec
    
    indx = np.where(Map == 1)
    pos_unsorted = np.array(indx)*disp_prec    
    pos_sorted = indx_sort(pos_unsorted,GMconfig)
    
    return pos_sorted

# Generate the disttable
def dist_table_gen(pos_sorted):
    # Record distances between all points
    reflectorNum = len(pos_sorted[0,:])
    dist_table = np.zeros([reflectorNum,reflectorNum]) # Full map
    
    for ii in range(reflectorNum):
        PtA = pos_sorted[:,ii]
        for jj in range(reflectorNum):
            PtB = pos_sorted[:,jj]
            disttmp = distmetrics(PtA,PtB)
            
            dist_table[jj,ii] = disttmp.copy()
            
    return dist_table
    
def ang_table_gen(pos_sorted):
    # Record distances between all points
    reflectorNum = len(pos_sorted[0,:])
    ang_table = np.zeros([reflectorNum,reflectorNum]) # Full map
    
    for ii in range(reflectorNum):
        PtA = pos_sorted[:,ii]
        for jj in range(reflectorNum):
            PtB = pos_sorted[:,jj]
            angtmp = angmetrics(PtA,PtB)
            
            ang_table[jj,ii] = angtmp
    
    return ang_table

# Generate Global Map based on reflector pos
def dispGlobalMap(reflector_pos,GMconfig,step=0):
    
    disp_prec = GMconfig.disp_prec
    simu_prec = GMconfig.simu_prec
    global_width = GMconfig.global_width
    global_length = GMconfig.global_length
    
    # Generate the global map
    GlobalMap = np.zeros([np.uint(global_width/disp_prec),np.uint(global_length/disp_prec)])
    
    reflector_pos_disp = np.array(np.round(reflector_pos/disp_prec),dtype='int')
    # Boundary handling
    indx = np.where(reflector_pos_disp[0,:] >= np.uint(global_width/disp_prec))
    reflector_pos_disp[0,indx] = np.uint(global_width/disp_prec)-1
    indx = np.where(reflector_pos_disp[1,:] >= np.uint(global_length/disp_prec))
    reflector_pos_disp[1,indx] = np.uint(global_length/disp_prec)-1
    
    _,reflectorNum = reflector_pos_disp.shape
    mask = 1
    for ii in range(reflectorNum):
        GlobalMap[reflector_pos_disp[0,ii],reflector_pos_disp[1,ii]] = mask
        mask += step
        
    return reflector_pos_disp,GlobalMap 
    # Notice that reflector_pos_disp is in a non-physical unit (pixels)
    
# Generate Map based oon reflector pos
def dispMap(reflector_pos,GMconfig,margin = 5,step=0):
    
    disp_prec = GMconfig.disp_prec
    simu_prec = GMconfig.simu_prec

    reflector_pos_disp = np.array(np.round(reflector_pos/disp_prec),dtype='int')
    
    width = np.round(np.nanmax(reflector_pos_disp[0,:]))+margin
    length = np.round(np.nanmax(reflector_pos_disp[1,:]))+margin
    width = np.uint(width)
    length = np.uint(length)
    
    # Generate the map
    Map = np.zeros([np.uint(width),np.uint(length)])
    
    _,reflectorNum = reflector_pos_disp.shape
    mask = 1
    for ii in range(reflectorNum):
        Map[reflector_pos_disp[0,ii],reflector_pos_disp[1,ii]] = mask
        mask += step
        
    return Map
    
# Generate Map based on reflector pos
# This function requires knowledge of the reflector pos
def dispMapW(reflector_pos,GMconfig,length,width):
    
    disp_prec = GMconfig.disp_prec
    simu_prec = GMconfig.simu_prec

    reflector_pos_disp = np.array(np.round(reflector_pos/disp_prec),dtype='int')
    
    # Generate the map
    Map = np.zeros([np.uint(width)+1,np.uint(length)+1])
    
    _,reflectorNum = reflector_pos_disp.shape
    for ii in range(reflectorNum):
        Map[reflector_pos_disp[0,ii],reflector_pos_disp[1,ii]] = 1
        
    return Map
