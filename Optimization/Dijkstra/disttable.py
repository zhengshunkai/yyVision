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
        ang = np.array(np.arctan2((PtA[0]-PtB[0]),(PtA[1]-PtB[1])),dtype='float')
        
    return ang

def angadjrange(ang):
    
    inflag = 0
    while inflag == 0:
        if ang > pi:
            ang = ang - 2*pi
        elif ang < -pi:
            ang = ang + 2*pi
        else:
            inflag = 1
    
    return ang
    
# Sort the control points
def indx_sort(pos_unsorted):
    
    # Sorting rule(s)
    # pos_unsorted_indx = pos_unsorted[0,:]*global_length+pos_unsorted[1,:] # line by line search
    # pos_unsorted_indx = np.sqrt(np.square(pos_unsorted[0,:])+np.square(pos_unsorted[1,:])) # Radius min to max search
    pos_unsorted_indx = pos_unsorted[0,:]*pos_unsorted[1,:]
    # pos_unsorted_indx = pos_unsorted[0,:]+pos_unsorted[1,:] # coordinates sum
   
    pos_ord = np.argsort(pos_unsorted_indx)
    pos_sorted = np.array(pos_unsorted[:,pos_ord],dtype='float')
    
    return pos_sorted
# Ideally, the distance between two points in the disttable shall represent the distance from the map (say points 1 and 100 shall be farther away from points 1 and 10)

# Generate position information from a map
def map2pos(Map):
    
    disp_prec = 1
    
    indx = np.where(Map == 1)
    pos_unsorted = np.array(indx)*disp_prec    
    pos_sorted = indx_sort(pos_unsorted)
    
    return pos_sorted

# Generate the disttable
def dist_table_gen(pos_sorted):
    # Record distances between all points
    CtrlPtsNum = len(pos_sorted[0,:])
    dist_table = np.zeros([CtrlPtsNum,CtrlPtsNum]) # Full map
    
    for ii in range(CtrlPtsNum):
        PtA = pos_sorted[:,ii]
        for jj in range(CtrlPtsNum):
            PtB = pos_sorted[:,jj]
            disttmp = distmetrics(PtA,PtB)
            
            dist_table[jj,ii] = disttmp.copy()
            
    return dist_table
    
def ang_table_gen(pos_sorted):
    # Record distances between all points
    CtrlPtsNum = len(pos_sorted[0,:])
    ang_table = np.zeros([CtrlPtsNum,CtrlPtsNum]) # Full map
    
    for ii in range(CtrlPtsNum):
        PtA = pos_sorted[:,ii]
        for jj in range(CtrlPtsNum):
            PtB = pos_sorted[:,jj]
            angtmp = angmetrics(PtA,PtB)
            angtmp2 = angmetrics(PtB,PtA)
            
            ang_table[jj,ii] = angtmp
    
    return ang_table
    
# Generate Map based oon CtrlPts pos
def dispMap(CtrlPts_pos,Width,Length,step=1):
    
    disp_prec = 1
    simu_prec = 1

    CtrlPts_pos_disp = np.array(np.round(CtrlPts_pos/disp_prec),dtype='int')
    
    # Generate the map
    Map = np.zeros([np.uint(Width),np.uint(Length)])
    
    _,CtrlPtsNum = CtrlPts_pos_disp.shape
    mask = 1
    for ii in range(CtrlPtsNum):
        Map[CtrlPts_pos_disp[0,ii],CtrlPts_pos_disp[1,ii]] = mask
        mask += step
        
    return Map