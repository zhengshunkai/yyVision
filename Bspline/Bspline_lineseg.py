## Bpline line segments approximation
# Author: Yang Yu
# Date: 2018/5/9
# Description: NA

import math
from math import pi
import numpy as np
import scipy as sp

from scipy import interpolate

import time

import B_spline as BS

import matplotlib.pyplot as plt


## Functions
class linsegment:
    Gb = []
    lineseg = []
    str_u = []
    end_u = []
    
    def __init__(self,Gb,lineseg,str_u,end_u):
        self.Gb = Gb
        self.lineseg = lineseg
        self.str_u = str_u
        self.end_u = end_u

def B_spline_lineseg_repr(n,k,Gb,simu_uprec=1e-2,tol=15e-2):
    
    simu_length = np.int(1/simu_uprec)+1
    
    us = np.array(range(simu_length))*simu_uprec
    M = us.shape[0]
    
    # In B-spline, knots (open and cardinal)
    knots = BS.knotsgen(n,k)
    
    # B-spline fitting (matrix)
    Approxline = BS.Bsplinefit(Gb,k,us)
    
    # Routine for line segment approximation
    str_indx = 0
    end_indx = simu_length-1
    
    line_segments = []
    
    # Figure
    # plt.figure()
    # plt.plot(Approxline[:,0],Approxline[:,1])
    
    Approxline = np.array(Approxline)
    while 1:
        # Line segment
        lineseg_length = end_indx - str_indx + 1
        str_pt = Approxline[str_indx,:]
        end_pt = Approxline[end_indx,:]
        str_u = us[str_indx]
        end_u = us[end_indx]
        
        # print(lineseg_length,str_indx,end_indx)
        
        delta_vector_lineseg = end_pt-str_pt
        Gb_lineseg = delta_vector_lineseg*np.array([[0,0],[1/3,1/3],[2/3,2/3],[1,1]])
        us_lineseg = us.copy()
        n_lineseg = len(Gb_lineseg)-1
        
        knots = BS.knotsgen(n_lineseg,k)
        lineseg = BS.Bsplinefit(Gb_lineseg,k,us)+str_pt
        lineseg = np.array(lineseg)
        
        # Find the max of minimum distance between a pt in lineseg and Approxline
        difflineval = np.zeros(lineseg_length)
        for ii in range(lineseg_length):
            diffline = np.sum(np.abs(Approxline[str_indx:end_indx+1,:] - lineseg[ii,:]),axis=1)
            difflineval[ii] = np.nanmin(diffline)
        
        difflinevalmax = np.max(difflineval)
        # print(difflinevalmax)
        
        if difflinevalmax > tol:
            lineseg_length = np.int(np.floor(lineseg_length/2))
            end_indx = str_indx + lineseg_length-1
        elif end_indx == simu_length-1:
            line_segment = linsegment(Gb_lineseg,lineseg,str_u,end_u)
            line_segments.append(line_segment)            
            break
        else:
            line_segment = linsegment(Gb_lineseg,lineseg,str_u,end_u)
            line_segments.append(line_segment)
            str_indx = end_indx
            end_indx = simu_length - 1
        
    return line_segments

def lineseg_disp(line_segments,Approxline):
    
    plt.figure()
    plt.plot(Approxline[:,0],Approxline[:,1])
    
    for ii in range(len(line_segments)):
        lineseg = line_segments[ii].lineseg
        plt.plot(lineseg[:,0],lineseg[:,1],'r')

def line_pt_dist(line_segments,pt):
    
    pt2linedist_min = np.inf
    intersect_pt_min = np.nan
    lineseg_min_indx = np.nan
    
    for indx in range(len(line_segments)):
        
        lineseg = line_segments[indx].lineseg
        pt2linedist,intersect_pt = lineseg_pt_dist(lineseg,pt)
        
        if pt2linedist < pt2linedist_min:
            pt2linedist_min = pt2linedist
            intersect_pt_min = intersect_pt.copy()
            lineseg_min_indx = indx
    
    return intersect_pt_min,lineseg_min_indx,pt2linedist_min

def lineseg_pt_dist(lineseg,pt):
    # Find the line corresponding to the lineseg
    # ax+by+c = 0: Av = 0 (A=[a,b,c],v=[x,y,1].T)
    A = np.zeros(3)
    Ap = np.zeros(3)
    if lineseg[0,0] == lineseg[-1,0]: # x constant
        A[0] = 0
        A[1] = 1
        A[2] = -lineseg[1,1]
    else: 
        A[0] = -(lineseg[-1,1]-lineseg[0,1])/(lineseg[-1,0]-lineseg[0,0])
        A[1] = 1
        A[2] = -(A[0]*lineseg[0,0]+A[1]*lineseg[0,1])
        
    # Find the distance from a point to the line
    # construct the perpendicular line
    Ap[0] = A[1]
    Ap[1] = -A[0]
    Ap[2] = -Ap[0]*pt[0]-Ap[1]*pt[1]
    
    # Intersection
    intersect_pt = line_line_intersect(A,Ap)
    pt2linedist = pt2pt_dist(intersect_pt,pt)
    
    offboundary = 0
    if lineseg[0,0] == lineseg[-1,0]: # x constant
        ydiff1 = intersect_pt[1] - lineseg[-1,1]
        ydiff0 = intersect_pt[1] - lineseg[0,1]
        if ydiff1*ydiff0 > 0: # All positive
            offboundary = 1
    else:
        xdiff1 = intersect_pt[0] - lineseg[-1,0]
        xdiff0 = intersect_pt[0] - lineseg[0,0]
        if xdiff1*xdiff0 > 0: # All positive
            offboundary = 1
    
    if offboundary == 1:
        pt_boundary_l = pt2pt_dist(pt,lineseg[0])
        pt_boundary_r = pt2pt_dist(pt,lineseg[1])
        if pt_boundary_l > pt_boundary_r:
            intersect_pt = lineseg[-1]
            pt2linedist = pt_boundary_r
        else:
            intersect_pt = lineseg[0]
            pt2linedist = pt_boundary_l
        
        intersect_pt = np.array(intersect_pt.T)
    
    return pt2linedist,intersect_pt
    
def line_line_intersect(A,Ap):
    
    AA = np.matrix([[A[0],A[1]],[Ap[0],Ap[1]]])
    bb = -np.matrix([A[2],Ap[2]])
    
    intersect_pt = np.array(np.linalg.pinv(AA)*bb.T)
    intersect_pt = intersect_pt.T[0]
    
    return intersect_pt

def pt2pt_dist(pt0,pt1):
    
    dist = np.sqrt(np.nansum(np.abs((pt0-pt1)**2)))
    
    return dist