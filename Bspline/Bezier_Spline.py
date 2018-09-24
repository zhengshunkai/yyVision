## Spline (Bezier line)
# Author: Yang Yu
# Date: 2018/4/23
# Description: Control model spline fitting (Bezier)

import math
from math import pi
import numpy as np
import scipy as sp
from scipy import misc as misc
import time
import matplotlib.pyplot as plt

## Functions
# Bezier matrix generation at u (n order)
def bezierMtK(n,u,k):
    
    Mtk = u**(k)*(1-u)**(n-k)*misc.comb(n,k)
    
    return Mtk

def bezierMb(us,n):
    
    length = len(us)
    Mb = np.matrix(np.zeros([length,n+1]))
    
    for k in range(n+1):
        for jj in range(length):
            u = us[jj]
            Mb[jj,k] = bezierMtK(n,u,k)
    
    return Mb
    
def bezierfit(Gb,Mb): # Cubic
    
    # Mb = bezierMb(us,n)
    Approline = Mb*Gb
    
    return Approline
    
# Bezier curve control points find rountine (pseudo inverse approach)
def lsqfit(Approxline,us,n=3):
    
    Mb = bezierMb(us,n)
    Mb_ = np.linalg.pinv(Mb)
    Gbe = Mb_ * Approxline
    
    return Gbe
    
## Cubic specific
def cubicbezierdMtK(u,k):
    
    if k == 0:
        dMtK = 3*u**2-6*u+3
    elif k == 1:
        dMtK = 9*u**2-12*u+3
    elif k == 2:
        dMtK = -9*u**2+6*u
    elif k == 3:
        dMtK = 3*u**2
    
    return dMtK

def cubicbezierdMb(us,n):
    
    length = len(us)
    dMb = np.matrix(np.zeros([length,n+1]))
    
    for k in range(n+1):
        for jj in range(length):
            u = us[jj]
            dMb[jj,k] = cubicbezierdMtK(u,k)
            
    return dMb
    
# Bezier curve derivative with respect to u
def Cubicbezierd(Gb,dMb): # Cubic only, numerical calculation
    
    Approlined = dMb*Gb
    
    return Approlined
    
def cubicbezierd2MtK(u,k):
    
    if k == 0:
        d2MtK = 6-6*u
    elif k == 1:
        d2MtK = -12+18*u
    elif k == 2:
        d2MtK = 6-18*u
    elif k == 3:
        d2MtK = 6*u
    
    return d2MtK

def cubicbezierd2Mb(us,n):
    
    length = len(us)
    d2Mb = np.matrix(np.zeros([length,n+1]))
    
    for k in range(n+1):
        for jj in range(length):
            u = us[jj]
            d2Mb[jj,k] = cubicbezierd2MtK(u,k)
            
    return d2Mb
    
def Cubicbezierd2(Gb,d2Mb): # Cubic only, numerical calculation
    
    Approlined2 = d2Mb*Gb
    
    return Approlined2

if __name__ == "__main__":
    ## Configuration
    simu_tprec = 1e-3
    simu_length_m = 1
    simu_length = np.int(simu_length_m/simu_tprec)
    n = 3 # Bezier order
    
    ## least square qbezier fit using penrose pseudoinverse
    V=np.array
    E,  W,  N,  S =  V((1,0)), V((-1,0)), V((0,1)), V((0,-1)) # Directions in global map
    
    # Find the control points
    # Area of simulation
    cw = 1
    ch = 1
    cps = V((0, 0)) # Start point
    cpe = V((cw, ch)) # End point, to the west
    
    # Gb = [cps,cps+cw*E/4,cps+cw*E/2,cpe] # Line
    Gb_delta = np.array([[0,0],[1/2,0],[1,1/2],[1,1]]) # delta percentage of cpe-cps vector
    delta_vector = cpe - cps
    Gb = Gb_delta*delta_vector # Curve control points (4 needed)
    
    us = V(range(simu_length+1))*simu_tprec
    
    Mb = bezierMb(us,n)
    Approline = bezierfit(Gb,Mb) #produces the approximation on the bezier curve at u in us. Fitting result
    
    CtrlPts_Est = lsqfit(Approline,us,n)
    # Catmull-Rom splines: always passing through the control points
    
    plt.figure()
    plt.plot(Approline[:,0],Approline[:,1])
    plt.plot(CtrlPts_Est[:,0],CtrlPts_Est[:,1],'--*')
    plt.grid()
    
    # Calculate the derivatives of the bezier curve (dx/dt,dy/dt)
    dMb = cubicbezierdMb(us,n)
    Approlinedt = Cubicbezierd(Gb,dMb)
    d2Mb = cubicbezierd2Mb(us,n)
    Approlinedt2 = Cubicbezierd2(Gb,d2Mb)
    
    plt.figure()
    plt.plot(Approlinedt[:,0],Approlinedt[:,1])
    plt.grid()
    
    plt.figure()
    plt.plot(Approlinedt2[:,0],Approlinedt2[:,1])
    plt.grid()
    
    # Mutliple points curve fitting (cubic,8 points fitting)
    Gb_delta_m = np.array([[0,0],[1/8,0],[1/4,1/8],[3/8,1/4],[1/2,3/8],[5/8,1/2],[3/4,5/8],[7/8,3/4],[1,7/8],[1,1]])
    Gb_m = Gb_delta_m*delta_vector # Curve control points (4 needed)
    
    n = len(Gb_delta_m)-1 # Requires a 9th order bezier to fit
    Mb = bezierMb(us,n)
    Approline_m = bezierfit(Gb_m,Mb)
    CtrlPts_Est_m = lsqfit(Approline_m,us,n)
    
    plt.figure()
    plt.plot(Approline_m[:,0],Approline_m[:,1])
    plt.plot(CtrlPts_Est_m[:,0],CtrlPts_Est_m[:,1],'--*')
    plt.grid()
    
    plt.show()