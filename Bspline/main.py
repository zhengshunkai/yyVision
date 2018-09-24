## Spline Pt distance
# Author: Yang Yu
# Date: 2018/4/27
# Description: Based on the assumption that rotation speed is constant (abstract)

import math
from math import pi
import numpy as np
import scipy as sp

from scipy import interpolate

import time

import B_spline as BS
import Bspline_lineseg as BSLS
import mindist as MD

import matplotlib.pyplot as plt

## Functions
#

if __name__ == "__main__":
    ## Configuration （B-spline）
    simu_uprec = 1e-3
    simu_xprec = simu_uprec
    simu_length = np.int(1/simu_uprec)
    k = 3 # B-spline (knots)
    
    E,  W,  N,  S =  np.array((1,0)), np.array((-1,0)), np.array((0,1)), np.array((0,-1)) # Directions in global map
    
    # Find the control points
    # Area of simulation
    cw = 15
    ch = 3
    cps = np.array((0, 0)) # Start point
    cpe = np.array((cw, ch)) # End point, to the west
    
    # Gb_delta = np.array([[0,0],[1/8,0],[1/4,1/8],[3/8,1/4],[1/2,3/8],[5/8,1/2],[3/4,5/8],[7/8,3/4],[1,7/8],[1,1]])
    Gb_delta = np.array([[0,0],[1/2,0],[1,1/2],[1,1]])
    delta_vector = cpe - cps
    Gb = Gb_delta*delta_vector # Curve control points (4 needed)
    n = len(Gb_delta)-1
    
    ## B-spline fitting
    us = np.array(range(simu_length+1))*simu_uprec
    M = us.shape[0]
    
    # In B-spline, knots (open and cardinal)
    knots = BS.knotsgen(n,k)
    
    # B-spline fitting (matrix)
    Mb = BS.BsplinefitMb(n,k,us,knots)
    Approxline = BS.Bsplinefit(Gb,k,us)
    
    # Find thex control points
    CtrlPts_Est = BS.Bspline_lsqfit(Approxline,us,n,k,knots)
    
    # # x step equalization
    xx,yy = BS.Approxline_xequalize(Approxline,simu_xprec)
    Approxline_xeq = np.zeros(Approxline.shape)
    Approxline_xeq[:,0] = xx
    Approxline_xeq[:,1] = yy
    
    plt.figure()
    plt.plot(Approxline[:,0],Approxline[:,1])
    plt.plot(CtrlPts_Est[:,0],CtrlPts_Est[:,1],'--*')
    plt.grid()
    
    ## B-spline line segment fitting
    tic = time.clock()
    line_segments = BSLS.B_spline_lineseg_repr(n,k,Gb)
    toc = time.clock()
    elapse3 = toc-tic
    
    pt = np.array([5,0.5]).T
    
    # lineseg = line_segments[1].lineseg
    # pt2linedist,intersect_pt = BSLS.lineseg_pt_dist(lineseg,pt)
    # perpendicular_lineseg = np.zeros([2,2])
    # perpendicular_lineseg[:,0] = np.array([pt[0],intersect_pt[0]])
    # perpendicular_lineseg[:,1] = np.array([pt[1],intersect_pt[1]])
    
    intersect_pt_min,lineseg_min_indx,pt2linedist_min = BSLS.line_pt_dist(line_segments,pt)
    perpendicular_lineseg = np.zeros([2,2])
    perpendicular_lineseg[:,0] = np.array([pt[0],intersect_pt_min[0]])
    perpendicular_lineseg[:,1] = np.array([pt[1],intersect_pt_min[1]])
    
    BSLS.lineseg_disp(line_segments,Approxline)
    plt.plot(perpendicular_lineseg[:,0],perpendicular_lineseg[:,1],'r--')
    plt.grid()
    
    # Examine the minization line
    Approxline_pt_diffline = np.array(np.abs(Approxline-pt))
    Approxline_pt_diff_val = np.sqrt(Approxline_pt_diffline[:,0]**2+Approxline_pt_diffline[:,1]**2)
    
    # plt.figure()
    # plt.plot(us,Approxline_pt_diff_val)
    
    u_init = 0.0
    itermax = 1000
    
    tic = time.clock()
    u = u_init
    for iter in range(itermax):
        # Newton's method
        u1,dist_grad,xy,dist = MD.update_u_newton(u,pt,n,k,Gb,delta_u = 1e-6,step = 1e-2)
        # print(u,dist_grad,dist,xy[1])
        u = u1
    toc = time.clock()
    
    print(u1,xy,dist)
    elapse0 = toc - tic
        
    # bracketing method
    tic = time.clock()
    min_u,min_xy,mindist = MD.mindist_bracketing(u_init,pt,n,k,Gb,delta_u = 1e-6,step = 1e-2,tol=1e-2)
    
    print(min_u,min_xy,mindist)
    
    toc = time.clock()
    elapse1 = toc - tic
    
    # Line u initialization
    tic = time.clock()
    u_init = MD.u_initialization(line_segments,pt)
    min_u,min_xy,mindist = MD.mindist_bracketing(u_init,pt,n,k,Gb,delta_u = 1e-6,step = 1e-2,tol=1e-2)
    toc = time.clock()
    elapse2 = toc - tic
    
    print(min_u,min_xy,mindist)
    print(elapse0,elapse1,elapse2)
    
    plt.show()