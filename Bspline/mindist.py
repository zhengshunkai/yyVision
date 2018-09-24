## Point to Bspline minimum distance cal.
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
import Bspline_lineseg as BSLS

## Functions
def u_initialization(line_segments,pt):
    
    intersect_pt_min,lineseg_min_indx,pt2linedist_min = BSLS.line_pt_dist(line_segments,pt)
    
    lineseg = line_segments[lineseg_min_indx].lineseg
    str_u = line_segments[lineseg_min_indx].str_u
    end_u = line_segments[lineseg_min_indx].end_u
    
    dist_str_intersectpt = BSLS.pt2pt_dist(lineseg[0],pt)
    dist_end_intersectpt = BSLS.pt2pt_dist(lineseg[-1],pt)
    
    u_est = (end_u-str_u)*(dist_str_intersectpt/(dist_end_intersectpt+dist_str_intersectpt))+str_u
    
    return u_est
    
def grad_cal(u,pt,n,d,Gb,delta_u = 1e-6): # Numerical approach
    
    us = np.array([u-delta_u,u,u+delta_u])
    indx = np.where(us>1)
    us[indx] = 1
    indx = np.where(us<0)
    us[indx] = 0 
    
    knots = BS.knotsgen(n,d)
    xy = BS.Bsplinefit(Gb,d,us)
    
    # Gradient
    # min((x-ptx)^2+(y-pty)^2)
    # grad u: 2(x-ptx)dx/du + 2(y-pty)dy/du
    delta_x = xy[-1,0]-xy[0,0]
    delta_y = xy[-1,1]-xy[0,1]
    delta_u = us[-1]-us[0]
    
    dist = np.sqrt((xy[1,0]-pt[0])**2+(xy[1,1]-pt[1])**2)
    dist_grad = (2*(xy[1,0]-pt[0])*delta_x/delta_u+2*(xy[1,1]-pt[1])*delta_y/delta_u)/(2*dist)
    
    xy = np.array(xy)[1]

    return dist_grad,dist,xy
    
def update_u_newton(u,pt,n,d,Gb,delta_u = 1e-6,step = np.nan):
    
    dist_grad,dist,xy = grad_cal(u,pt,n,d,Gb,delta_u)
    
    if np.isnan(step):
        u1 = u - dist/dist_grad # Newton's method tends to be instable
        if u1 < 0:
            u1 = 0
        elif u1 > 1:
            u1 = 1
    else: # Fixed step
        u1 = u - np.sign(dist_grad)*step
    
    return u1,dist_grad,xy,dist
# The above method yields a lot to be desired. First, boundaries u=0 and u=1 might be surpassed by a delta. Second, accuracy limited in both Newton's method and the fixed step method (which also takes longer to converge)

# A bracketing approach
def mindist_bracketing(u,pt,n,d,Gb,delta_u = 1e-6,step = 1e-2,tol = 1e-3):
    
    bracket = np.ones(2)*np.nan
    bracket_xy = np.ones([2,2])*np.nan
    
    # Initialization: finding the first bracket via Newton/fixed step method
    dist_grad,dist,xy = grad_cal(u,pt,n,d,Gb,delta_u)
    
    if dist_grad >= 0: # the zero case is trivia, but for completeness' sake
        bracket[1] = u
        bracket_xy[1] = xy
        dir = 1
    elif dist_grad < 0:
        bracket[0] = u
        bracket_xy[0] = xy
        dir = 0
    
    u = u - np.sign(dist_grad)*step # first iteration
    while 1: # Start searching the other boundary of the bracket
        
        dist_grad,dist,xy = grad_cal(u,pt,n,d,Gb,delta_u)
        
        if dist_grad < 0:
            bracket[0] = u # Update
            bracket_xy[0] = xy
            if dir == 1:
                break
            else:
                u = u - np.sign(dist_grad)*step
                
        elif dist_grad >= 0:
            bracket[1] = u
            bracket_xy[1] = xy
            if dir == 0:
                break
            else:
                u = u - np.sign(dist_grad)*step
    
    # Examine the center of the bracket (iteration)
    while 1:
        
        dist_r = BSLS.pt2pt_dist(bracket_xy[1],pt)
        dist_l = BSLS.pt2pt_dist(bracket_xy[0],pt)
        
        u_center = (bracket[1]+bracket[0])/2
        dist_grad,dist,xy = grad_cal(u_center,pt,n,d,Gb,delta_u)
        
        if np.abs(dist_l-dist_r)<=tol:
            min_u = u_center
            mindist = dist
            min_xy = xy
            break
        else:
            if dist_grad < 0:
                bracket[0] = u_center
                bracket_xy[0] = xy
            elif dist_grad >= 0:
                bracket[1] = u_center
                bracket_xy[1] = xy
        
        # print(bracket[0],bracket_xy[0],dist_l)
        # print(bracket[1],bracket_xy[1],dist_r)
        
    return min_u,min_xy,mindist
# The above method yields a lot to be desired. First, boundaries u=0 and u=1 might be surpassed by a delta. Second, accuracy limited in both Newton's method and the fixed step method (which also takes longer to converge)