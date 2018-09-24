## Spline (B-spline)
# Author: Yang Yu
# Date: 2018/4/23
# Description: Control model spline fitting (B-spline)

import math
from math import pi
import numpy as np
import scipy as sp
from scipy import misc as misc
import time

from deline import deline2, diff1D

import matplotlib.pyplot as plt

## Bspline generation routine
class SplinePath:
    Gb = []
    k = []
    
    def __init__(self,Gb=[],k=[]):
        self.Gb = Gb
        self.k = k
        
    def copy(self):
        SP_copy = SplinePath()
        SP_copy.Gb = self.Gb.copy()
        SP_copy.k = self.k
        
        return SP_copy

# Knots generation
def knotsgen(n,k):
    # Open and cardinal
    p = n+k+1
    knots = np.zeros(p)
    indx = 0
    for ii in range(p):
        if ii < k:
            knots[ii] = indx
        elif ii > p-k-1:
            knots[ii] = indx+1
        else:
            indx = indx+1
            knots[ii] = indx
    
    knots = knots/(indx+1)
    
    return knots

# knots generation with alpha
def knotsgen_alpha(n,k,alpha=0.25):
    # Open and cardinal
    p = n+k+1
    knots = np.zeros(p)
    indx = 0
    indxmax = p-2*k+1
    
    for ii in range(p):
        if ii < k:
            knots[ii] = indx
        elif ii > p-k-1:
            knots[ii] = indxmax
        else:
            indx = indx+1
            knots[ii] = indx
    
    knots = knots/(indx+1)
    
    for ii in range(p):
        if ii < k:
            pass
        elif ii > p-k-1:
            pass
        else:
            if ii<p/2:
                knots[ii] = knots[ii]*alpha
            else:
                knots[ii] = 1-(1-knots[ii])*alpha
    
    return knots
    
def knotsgen_beta(n,k,beta=0.5):
    # Open and cardinal
    p = n+k+1
    knots = np.zeros(p)
    indx = 0
    indxmax = p-2*k+1
    
    for ii in range(p):
        if ii < k:
            knots[ii] = indx
        elif ii > p-k-1:
            knots[ii] = indxmax
        else:
            indx = indx+1
            knots[ii] = indx
    
    knots = knots/(indx+1)
    alpha = beta/knots[k+1]
    for ii in range(p):
        if ii < k:
            pass
        elif ii > p-k-1:
            pass
        else:
            if ii<p/2:
                knots[ii] = knots[ii]*alpha
            else:
                knots[ii] = 1-(1-knots[ii])*alpha
    
    return knots
    
def knotsgen2(n,k,r1=0.5):
    # Open and cardinal
    p = n+k+1
    knots = np.zeros(p)
    indx = 0
    
    for ii in range(p):
        if ii < k:
            knots[ii] = indx
        elif ii > p-k-1:
            knots[ii] = indx+r1
        elif ii == k:
            indx = indx+r1
            knots[ii] = indx
        else:
            indx = indx+1
            knots[ii] = indx
            
    knots = knots/knots[-1]

    return knots
    
# Basic functions generations (require knots function)
# Ref: http://web.mit.edu/hyperbook/Patrikalakis-Maekawa-Cho/node16.html
def BsplineBnk(i,k,u,knots):
    
    u_delta = 1e-10
    
    ui = knots[i]
    ui1 = knots[i+1]
    
    if k == 1:
        # Boundary control
        if u == 1: # The end 
            u = u-u_delta
            
        if u>=ui and u<ui1:
            Bik = 1
        else:
            Bik = 0
    else:
        Bik_1 = BsplineBnk(i,k-1,u,knots)
        Bi1k_1 = BsplineBnk(i+1,k-1,u,knots)
        uik = knots[i+k]
        uik_1 = knots[i+k-1]
        
        if uik_1 == ui:
            alpha = 0
        else:
            alpha = (u-ui)/(uik_1-ui)
        
        if uik == ui1:
            beta = 0
        else:
            beta = (uik-u)/(uik-ui1)
        
        Bik = alpha*Bik_1 + beta*Bi1k_1
        
    return Bik

# Derivative
def BsplineBnk_derivative(i,k,u,knots):
    
    u_delta = 1e-10
    
    ui = knots[i]
    ui1 = knots[i+1]
    
    if u == 1:
        u = u - u_delta
    
    # Calculate the splines
    Bik_1 = BsplineBnk(i,k-1,u,knots)
    Bi1k_1 = BsplineBnk(i+1,k-1,u,knots)
        
    uik = knots[i+k]
    uik_1 = knots[i+k-1]
        
    if uik_1 == ui:
        alpha = 0
    else:
        alpha = (k-1)/(uik_1-ui)
        
    if uik == ui1:
        beta = 0
    else:
        beta = -(k-1)/(uik-ui1)
        
    deBik = alpha*Bik_1 + beta*Bi1k_1
        
    return deBik
    
def BsplineBnk_derivative2(i,k,u,knots):
    
    u_delta = 1e-10
    
    ui = knots[i]
    ui1 = knots[i+1]
    
    if u == 1:
        u = u - u_delta
    
    # Calculate the splines
    # 2nd order
    Bik_1 = BsplineBnk_derivative(i,k-1,u,knots)
    Bi1k_1 = BsplineBnk_derivative(i+1,k-1,u,knots)
        
    uik = knots[i+k]
    uik_1 = knots[i+k-1]
        
    if uik_1 == ui:
        alpha = 0
    else:
        alpha = (k-1)/(uik_1-ui)
        
    if uik == ui1:
        beta = 0
    else:
        beta = -(k-1)/(uik-ui1)
        
    de2Bik = alpha*Bik_1 + beta*Bi1k_1
        
    return de2Bik
    
def BsplinefitMb(n,k,us,knots):
    
    Mb = np.matrix(np.zeros([len(us),n+1]))
    
    for ii in range(len(us)):
        u = us[ii]
        for nn in range(n+1):
            Mb[ii,nn] = BsplineBnk(nn,k,u,knots)
            
    return Mb
    
def BsplinefitMb_derivative(n,k,us,knots):
    
    deMb = np.matrix(np.zeros([len(us),n+1]))
    
    for ii in range(len(us)):
        u = us[ii]
        for nn in range(n+1):
            deMb[ii,nn] = BsplineBnk_derivative(nn,k,u,knots)
            
    return deMb
    
def BsplinefitMb_derivative2(n,k,us,knots):
    
    de2Mb = np.matrix(np.zeros([len(us),n+1]))
    
    for ii in range(len(us)):
        u = us[ii]
        for nn in range(n+1):
            de2Mb[ii,nn] = BsplineBnk_derivative2(nn,k,u,knots)
            
    return de2Mb

# Fit spline
def Bsplinefit(Gb,k,us):
    
    n = len(Gb) - 1
    knots = knotsgen(n,k)
    
    Mb = BsplinefitMb(n,k,us,knots)
    Approline = Mb*Gb
    
    return Approline
    
# Fit spline derivative
def Bsplinefit_derivative(Gb,k,us):
    
    n = len(Gb) - 1
    knots = knotsgen(n,k)
    
    deMb = BsplinefitMb_derivative(n,k,us,knots)
    deApproline = deMb*Gb
    
    return deApproline
    
# Fit spline derivative
def Bsplinefit_derivative2(Gb,k,us):
    
    n = len(Gb) - 1
    knots = knotsgen(n,k)
    
    de2Mb = BsplinefitMb_derivative2(n,k,us,knots)
    de2Approline = de2Mb*Gb
    
    return de2Approline
    
# B-spline curve control points finding rountine (pseudo inverse approach)
def Bspline_lsqfit(Approxline,us,n,d,knots):
    
    Mb = BsplinefitMb(n,d,us,knots)
    Mb_ = np.linalg.pinv(Mb)
    Gbe = Mb_ * Approxline
    
    return Gbe

def Approxline_xequalize(Approxline,simu_xprec):
    x = np.array(Approxline)[:,0]
    y = np.array(Approxline)[:,1]
    
    cw = x[-1]-x[0]
    xx = np.arange(np.int(1/simu_xprec)+1)*simu_xprec*cw
    
    f = sp.interpolate.interp1d(x,y)
    
    yy = f(xx)
    
    return xx,yy

def curvature_cal_numerical(trajectory,mode=0):
    
    x = vecMat2Array(trajectory[:,0])
    y = vecMat2Array(trajectory[:,1])
    
    dx = diff1D(x)
    dy = diff1D(y)
    ddx = diff1D(dx)
    ddy = diff1D(dy)
    
    if mode == 1:
        curvature = np.abs(dx*ddy-dy*ddx)/(dx**2+dy**2)**(3/2)
    else:
        curvature = (dx*ddy-dy*ddx)/(dx**2+dy**2)**(3/2)
    
    curvature[-1] = curvature[-3]
    curvature[-2] = curvature[-3]
    
    return curvature
    
def curvature_cal(Gb,k,us,mode=0):
    
    deApproline = Bsplinefit_derivative(Gb,k,us)
    de2Approline = Bsplinefit_derivative2(Gb,k,us)
    
    dx = vecMat2Array(deApproline[:,0])
    dy = vecMat2Array(deApproline[:,1])
    ddx = vecMat2Array(de2Approline[:,0])
    ddy = vecMat2Array(de2Approline[:,1])
    
    if mode == 1:
        curvature = np.abs(dx*ddy-dy*ddx)/(dx**2+dy**2)**(3/2)
    else:
        curvature = (dx*ddy-dy*ddx)/(dx**2+dy**2)**(3/2)
    
    return curvature
    
def vecMat2Array(vecM): # vector from matrix format to array
    vecA = np.array(vecM)
    vecA = np.reshape(vecA,[1,-1])[0]
    
    return vecA

def spline_eqdist_us(Gb,k,dist_avg=1e-2): # Making the distance around 1e-2 m between points
    
    # Calculate the length of the spline by measuring the distances between Gb
    pathdist = 0
    for ii in range(len(Gb)-1):
        pt1 = Gb[ii]
        pt2 = Gb[ii+1]
        pathdist = pathdist + np.sqrt(np.sum((pt1-pt2)**2))
        
    us_prec = dist_avg/pathdist/10
    
    str_u = 0
    u = np.array([str_u])
    
    tol = dist_avg/1e2
    endflag = 0
    us = []
    us_dist = []
    us.append(u)
    while endflag == 0:
        pt1 = Bsplinefit(Gb,k,u)
        pt1 = np.array(pt1)
        
        distflag = 0
        while distflag == 0:
            
            u = u+us_prec
            
            if u >= 1:
                u = np.array([1])
                pt2 = Bsplinefit(Gb,k,u)
                pt2 = np.array(pt2)
                dist = np.sqrt(np.sum((pt1-pt2)**2))
                
                us_dist.append(dist)
                us.append(u)

                endflag = 1
                break
                
            pt2 = Bsplinefit(Gb,k,u)
            pt2 = np.array(pt2)
            
            dist = np.sqrt(np.sum((pt1-pt2)**2))
            if (dist-dist_avg)<tol:
                distflag = 0
            else:
                us_dist.append(dist)
                us.append(u)
                distflag = 1
    
    if us_dist[-1] < dist_avg/10:
        L = len(us)
        us[-2] = us[-1]
        us = us[0:L-1]
        us_dist[-2] = us_dist[-2] + us_dist[-1]
        us_dist = us_dist[0:L-1]
        
    us = np.reshape(np.array(us),[1,-1])[0]
    us_dist = np.array(us_dist)
    
    return us,us_dist

def spline_curvature_display(SPpath,simu_prec):
    
    simu_length = np.int(1/simu_prec)
    us = np.array(range(simu_length+1))*simu_uprec
    
    # Fit spline
    Gb = SP_path.Gb
    k = SP_path.k
    path = Bsplinefit(Gb,k,us)
    
    # Spline derivative
    # path_derivative = Bsplinefit_derivative(Gb,k,us)
    # path_derivative2 = Bsplinefit_derivative2(Gb,k,us)
    
    path_curvature = curvature_cal(Gb,k,us,mode=1)
    de_path_curvature = deline2(us,path_curvature)
    
    # Plot
    plt.figure()
    plt.plot(path[:,0],path[:,1])
    plt.plot(Gb[:,0],Gb[:,1],'r--')

    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(us,path_curvature)
        
    plt.subplot(2,1,2)
    plt.plot(us,de_path_curvature)

## Routine
if __name__ == "__main__":
    
    simu_uprec = 1e-3
    simu_xprec = simu_uprec
    simu_length = np.int(1/simu_uprec)
    k = 4 # B-spline (knots)
    
    V=np.array
    E,  W,  N,  S =  V((1,0)), V((-1,0)), V((0,1)), V((0,-1)) # Directions in global map
    
    # Find the control points
    # Area of simulation
    cw = 1
    ch = 1
    cps = V((0, 0)) # Start point
    cpe = V((cw, ch)) # End point, to the west
    
    # Gb_delta = np.array([[0,0],[1/8,0],[1/4,1/8],[3/8,1/4],[1/2,3/8],[5/8,1/2],[3/4,5/8],[7/8,3/4],[1,7/8],[1,1]])
    # Gb_delta = np.array([[0,0],[1/4,0],[1,3/4],[1,1],[1,1],[1,1],[1,1]])
    Gb_delta = np.array([[0,0],[1000,0],[1666.667,666.667],[2000,2000],[2000,4000]])/1000
    # Gb_delta = np.array([[0,0],[1/4,0],[[3/4,1],[1,1]]) # S shape
    delta_vector = cpe - cps
    Gb = Gb_delta*delta_vector # Curve control points (4 needed)
    n = len(Gb_delta)-1
    
    us = V(range(simu_length+1))*simu_uprec
    
    # knots generations
    knots = knotsgen(n,k)
    
    # General function generation
    Mb = BsplinefitMb(n,k,us,knots)
    deMb = BsplinefitMb_derivative(n,k,us,knots)
    de2Mb = BsplinefitMb_derivative2(n,k,us,knots)
    
    # plt.figure()
    # for nn in range(n+1):
    #     plt.plot(us,Mb[:,nn])
    #     
    # plt.figure()
    # for nn in range(n+1):
    #     plt.plot(us,deMb[:,nn])
    #     
    # plt.figure()
    # for nn in range(n+1):
    #     plt.plot(us,de2Mb[:,nn])
        
    # Fit spline
    Approline = Bsplinefit(Gb,k,us)
    
    # Spline derivative
    deApproline = Bsplinefit_derivative(Gb,k,us)
    de2Approline = Bsplinefit_derivative2(Gb,k,us)
    
    # Derivative: numerical simulation
    dxdu = deline2(us,Approline[:,0])
    dydu = deline2(us,Approline[:,1])
    dx2du2 = deline2(us,dxdu)
    dy2du2 = deline2(us,dydu)
    
    plt.figure()
    plt.plot(Approline[:,0],Approline[:,1])
    plt.plot(Gb[:,0],Gb[:,1],'--*')
    
    plt.figure()
    plt.plot(us,deApproline[:,0])
    plt.plot(us,dxdu,'r--')
    plt.figure()
    plt.plot(us,deApproline[:,1])
    plt.plot(us,dydu,'r--')
    
    plt.figure()
    plt.plot(us,de2Approline[:,0])
    plt.plot(us,dx2du2,'r--')
    plt.figure()
    plt.plot(us,de2Approline[:,1])
    plt.plot(us,dy2du2,'r--') # Notice the abnormality at the end of the spline, which is due to edge handling in deline2.
    
    # Curvature calculation
    curvature_numerical = curvature_cal_numerical(Approline)
    curvature_analytical = curvature_cal(Gb,k,us)
    
    plt.figure()
    plt.plot(us,curvature_numerical)
    plt.plot(us,curvature_analytical,'r--')
    
    plt.show()