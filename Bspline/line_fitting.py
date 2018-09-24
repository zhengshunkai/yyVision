## Line fitting
# Author: Yang Yu
# Date: 2018/5/11
# Description: NA

import math
from math import pi
import numpy as np

import time

import matplotlib.pyplot as plt

import meshop as MP
import centerop as CP

## Functions
V = np.array

class lineseg:
    # Line defined as ax+by+1 = 0
    a = []
    b = []
    c = []
    cps = [] # From 0 to 1
    cpe = []
    
    def __init__(self,a,b,c,cps,cpe):
        self.a = a
        self.b = b
        self.c = c
        self.cps = cps
        self.cpe = cpe
        
class line:
    a = []
    b = []
    c = []
    
    def __init__(self,a=0,b=0,c=0):
        self.a = a
        self.b = b
        self.c = c

def line2lineseg(l,cps,cpe):
    
    lseg = lineseg(l.a,l.b,l.c,cps,cpe)
    
    return lseg
    
def linesegconst(pt1,pt2):
    
    l = lineconst(pt1,pt2)
    lseg = line2lineseg(l,pt1,pt2)
    
    return lseg
    
def plot_linseg(lseg):
    
    cps = lseg.cps
    cpe = lseg.cpe
    
    x = np.array([cps[0],cpe[0]])
    y = np.array([cps[1],cpe[1]])
    
    plt.plot(x,y)
    
def lineseg2line(lseg):
    
    l = line(lseg.a,lseg.b,lseg.c)
    
    return l

def linelinesegintersect(lseg,l):
    
    l1 = lineseg2line(lseg)
    pti = lineintersect(l1,l)
    
    if len(pti) == 0: # Parallel case
        if ptonline(lseg.cps,l) == 1: # on line
            pti = np.array([lseg.cps,lseg.cpe]) # Overlap
        else:
            pti = []
        return pti
        
    flag = ptonlseg(pti,lseg)
    if flag == 0: # no intersection
        pti = np.array([])
        return pti
    elif flag == 1:
        return pti

def lseglsegintersect(lseg1,lseg2):
    
    l1 = lineseg2line(lseg1)
    pti = linelinesegintersect(lseg2,l1) # len = 2 or 1
    
    if len(pti) == 0: # Empty
        return pti
    elif len(pti) == 4: # Overlap
        flag1 = ptonlseg(lseg1.cps,lseg2)
        flag2 = ptonlseg(lseg1.cpe,lseg2)
        flag3 = ptonlseg(lseg2.cps,lseg1)
        flag4 = ptonlseg(lseg2.cpe,lseg1)
        
        if (flag1 == 1) and (flag2 != 1):
            pti = np.array([lseg1.cps,lseg2.cpe])
        elif (flag2 == 1) and (flag1 != 1):
            pti = np.array([lseg2.cps,lseg1.cpe])
        elif (flag1 == 1) and (flag2 == 1):
            pti = np.array([lseg1.cps,lseg1.cpe])
        elif (flag3 == 1) and (flag4 == 1):
            pti = np.array([lseg2.cps,lseg2.cpe]) # len = 4
        else:
            pti = [] # len = 0
    else:
        flag0 = ptonlseg(pti,lseg1) # Check if the point is also on lseg1
        if flag0 != 1: 
            pti = []
    
    return pti
    
def ptlinedist(pt,l):
    # Find the perpendicular line from pt to l
    lp = pendilineconst(l,pt)
    
    # Find the intersect
    pti = lineintersect(l,lp)
    
    # Find the distance between pt and pti
    dist = np.sqrt(np.sum((pt-pti)**2))
    
    return dist,pti

def pt2linesegdist(pt,lseg):
    
    l = lineseg2line(lseg)
    dpl,pti = ptlinedist(pt,l)
    
    flag = ptonlseg(pti,lseg)
    if flag == 0: # Not on the line segment
        d1 = pt2dist(pti,cps)
        d2 = pt2dist(pti,cpe)
        if d1 >= d2:
            d_min = d2
            pti = cps.copy()
        elif d2 > d1:
            d_min = d1
            pti = cpe.copy()
            
        return d_min,pti
    elif flag == 1: # On the line segment
        
        return dpl,pti
    
def ptonline(pt,l):
    
    flag = 1
    tol = 1e-7
    
    dpl,pti = ptlinedist(pt,l)
    if dpl >= tol:
        flag = 0
    
    return flag
    
def ptonlseg(pt,lseg):
    
    tol = 1e-7
    
    l = lineseg2line(lseg)
    dpl,pti = ptlinedist(pt,l)
    if dpl >= tol:
        flag = 0
    
    x = pt[0]
    y = pt[1]
    
    xmin = np.min([lseg.cps[0],lseg.cpe[0]])
    xmax = np.max([lseg.cps[0],lseg.cpe[0]])
    ymin = np.min([lseg.cps[1],lseg.cpe[1]])
    ymax = np.max([lseg.cps[1],lseg.cpe[1]])
    
    if np.abs(xmin-xmax) < tol:
        if (y < ymin) or (y > ymax):
            flag = 0
        else:
            flag = 1
    else:
        if (x < xmin) or (x > xmax):
            flag = 0
        else:
            flag = 1
    
    return flag
    
def lineconst(pt1,pt2):
    
    tol = 1e-5
    if np.abs(pt1[0]-pt2[0])<tol: # pt1 x == pt2 x
        a = V(1.0)
        b = V(0.0)
        c = -pt1[0]
    else:
        b = V(1.0)
        a = -(pt2[1]-pt1[1])/(pt2[0]-pt1[0])
        c = -a*pt1[0]-b*pt1[1]
        
    line_C = line(a,b,c)
    
    return line_C
    
def lineconst_k(pt1,angle,radflag = 0):
    
    if radflag == 0:
        angle = angle/180*pi
    
    pt2 = pt1 + np.array([np.cos(angle),np.sin(angle)])
    
    line_C = lineconst(pt1,pt2)
    
    return line_C

def pendilineconst(Lin,pt):
    Lp = line(0,0,0)
    if Lin.a == 0 and Lin.b == 0:
        return Lp
    elif Lin.a == 0:
        Lp.a = V(1.0)
        Lp.b = V(0.0)
        Lp.c = -pt[0]
    elif Lin.b == 0:
        Lp.a = V(0.0)
        Lp.b = V(1.0)
        Lp.c = -pt[1]
    else:
        Lp.a = -Lin.b/Lin.a
        Lp.b = V(1.0)
        Lp.c = -(Lp.a*pt[0]+Lp.b*pt[1])
        
    return Lp

def linedirvec(l):
    dirvec = np.array([l.b,-l.a])
    if l.b < 0:
        dirvec = -dirvec
    elif l.b == 0:
        if -l.a < 0:
            dirvec = -dirvec
        elif l.a == 0:
            print('Warning: a=b=0.')
            return -1
    
    dirvec = dirvec/veclength(dirvec)
    
    return dirvec

def veclength(v):
    
    vlen = np.sqrt(np.sum(v[0]**2+v[1]**2))
    
    return vlen
    
def vecnorm(v):
    
    vlen = veclength(v)
    v = v/vlen
    
    return v

def lineintersect(L1,L2):
    
    intersect = []
    
    if L1.a/L1.b == L2.a/L2.b: # Parallel
        return intersect
    else:
        A = np.matrix([[L1.a,L1.b],[L2.a,L2.b]])
        B = -np.matrix([[L1.c,L2.c]]).T
        
        v = np.linalg.pinv(A)*B
        v = np.array(v.T)[0]
    
        return v
        
def ptsdist(pt1,pt2):
    
    dist = np.sqrt(np.sum((pt1-pt2)**2))
    
    return dist

def dispLineseg(lS,simu_prec= 1e-3):
    
    tol = simu_prec/2
    
    a = lS.a
    b = lS.b
    c = lS.c
    
    Vord = np.int(1/simu_prec)
    Hord = np.int(1/simu_prec)
    # C = CP.center(Hord,Vord)
    
    mg = MP.meshgridMN(Vord,Hord,simu_prec,simu_prec)
    
    # Filtering
    if len(lS.edgex) > 0:
        if lS.edgex[0] != lS.edgex[1]:
            indx = np.where((np.abs(mg.meshx*a+mg.meshy*b+c)<tol)*(mg.meshx>=lS.edgex[0])*(mg.meshx<=lS.edgex[1]))
            x = mg.meshx[indx]
            y = mg.meshy[indx]
        else:
            indx = np.where((np.abs(mg.meshx*a+mg.meshy*b+c)<tol)*(mg.meshy>=lS.edgey[0])*(mg.meshy<=lS.edgey[1]))
            x = mg.meshx[indx]
            y = mg.meshy[indx]
    
    return x,y

def observe_lineseg(obs_pos,lS,std_dist,std_pendi,observation_num):
    
    # Observed ground truth
    x,y = dispLineseg(lS,simu_prec=1e-2)
    x_obs_GT = np.zeros(observation_num)
    y_obs_GT = np.zeros(observation_num)
    for ii in range(observation_num):
        x_obs_GT[ii] = x[0]+(x[-1]-x[0])/(observation_num-1)*ii
        y_obs_GT[ii] = y[0]+(y[-1]-y[0])/(observation_num-1)*ii
    
    # For convenience
    ptx = obs_pos[0]
    pty = obs_pos[1]
    
    x_obs = np.zeros(observation_num)
    y_obs = np.zeros(observation_num)
    for ii in range(observation_num):
        
        # Calculate the angle
        ang = np.arctan((pty-y_obs_GT[ii])/(ptx-x_obs_GT[ii]))/pi*180
        delta_dist = np.random.randn(1)*std_dist
        delta_pendi = np.random.randn(1)*std_pendi
        
        delta_x,delta_y = RotXY(delta_dist,delta_pendi,ang)
        # x_obs[ii] = x[0]+delta_x
        # y_obs[ii] = y[0]+delta_y
        x_obs[ii] = x_obs_GT[ii]+delta_x
        y_obs[ii] = y_obs_GT[ii]+delta_y
    
    return x_obs,y_obs

def line_fitting(x,y):

    if np.nanstd(x) < 1e-2:
        b = 0.0 # Form x = h: a=1,b=0,c=-h
        a = 1.0
        c = -(np.nanmean(y))
    else:
        # General form: ax+by+c = 0
        # y = kx+h: b=1,k=a,h=-c
        # Construction
        b = 1.0
        X = np.matrix([x,np.ones(len(x))]).T
        Y = np.matrix([y]).T
        A = np.linalg.pinv(X.T*X)*X.T*Y
        A = np.array(A.T)
        a = -A[0][0]
        c = -A[0][1]
    
    line_f = line(a,b,c)
    
    return line_f

def inprod(v):
    
    v = np.matrix(v)
    v_mean = np.sqrt(np.array(v*v.T)[0][0])
    
    return v_mean

def lineseg_measurement(x_obs,y_obs,obs_pos,lS1_length,x=[],y=[],pflag=1,debug=0):
    
    # Given observer position, and obervationss
    # Find line segment
    
    # Filter each x_obs column
    # Simple mean operation (Consider RANSAC)
    M,N = x_obs.shape
    x_obs_sel = np.zeros(N)
    y_obs_sel = np.zeros(N)
    for ii in range(N):
        x_obs_col = x_obs[:,ii]
        y_obs_col = y_obs[:,ii]
        x_obs_sel[ii] = np.nanmean(x_obs_col)
        y_obs_sel[ii] = np.nanmean(y_obs_col)
        
    # Fit line
    line_f = line_fitting(x_obs_sel,y_obs_sel)
    
    # Project obseved points to line_f
    proj_line_x = np.zeros(len(x_obs_sel))
    proj_line_y = np.zeros(len(y_obs_sel))
    
    for ii in range(len(x_obs_sel)):
        pt1 = obs_pos
        pt2 = np.array([x_obs_sel[ii],y_obs_sel[ii]])
        
        if pflag == 1:
            line_obs = lineconst(pt1,pt2)
            intersect_x,intersect_y = lineintersect(line_obs,line_f)
        elif pflag == 0:
            line_f_p = pendilineconst(line_f,pt2)
            intersect_x,intersect_y = lineintersect(line_f_p,line_f)
        
        proj_line_x[ii] = intersect_x
        proj_line_y[ii] = intersect_y
        
    if debug == 1:
        plt.figure()
        if len(x) > 0:
            plt.plot(x,y)
        plt.plot(proj_line_x,proj_line_y,'r--')
        plt.grid()
    
    # Fit in the line segment
    lineseg_center = np.zeros(2)
    lineseg_center[0] = (proj_line_x[-1]+proj_line_x[0])/2
    lineseg_center[1] = (proj_line_y[-1]+proj_line_y[0])/2
    
    # Edge calculation
    if proj_line_x[-1] != proj_line_x[0]:
        line_f_vector = np.array([proj_line_x[-1]-proj_line_x[0],proj_line_y[-1]-proj_line_y[0]])
        line_f_vector = line_f_vector/inprod(line_f_vector)
    else:
        line_f_vector = np.array([0,1])
    
    lineseg_left_edge = -line_f_vector*lS1_length/2+lineseg_center
    lineseg_right_edge = line_f_vector*lS1_length/2+lineseg_center
    
    return lineseg_center,lineseg_left_edge,lineseg_right_edge

## Affine model rotation
def AffineTransform(pos1,AffineMat):
    
    pos1 = np.reshape(pos1,[2,1])
    
    pos1e = vec_extend(pos1)
    pos1eR = AffineMat*np.matrix(pos1e)
    
    pos1R = np.zeros([2,1])
    pos1R[0,:] = pos1eR[0,:]
    pos1R[1,:] = pos1eR[1,:]
    pos1R = np.reshape(np.array(pos1R),[1,-1])[0]
    
    return pos1R
    
def vec_extend(v):
    
    M,N = v.shape
    ve = np.ones([3,N])
    ve[0,:] = v[0,:]
    ve[1,:] = v[1,:]
    
    return ve
    
def affineRotA(theta,radflag = 0, extendflag = 1): # rad
    if radflag == 0:
        theta = theta/180*pi
    
    if extendflag == 0:
        RotA = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    else:
        RotA = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    
    return RotA
    
def ptRot(pt,angle):
    
    A = affineRotA(angle,radflag = 0, extendflag = 1)
    ptR = AffineTransform(pt,A)
    
    return ptR
    
def tangent_vec(v):
    
    v = vecnorm(v)
    ang = np.arctan2(v[1],v[0])/pi*180
    
    ang1 = ang + 90
    
    v1 = np.array([np.cos(ang1/180*pi),np.sin(ang1/180*pi)])
    
    return v1

def RotMatrix(theta,degflag = 1): 
    
    theta = theta/180*pi # Theta in degree/rad
    A = np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    
    return A
    
def RotXY(x,y,theta):
    
    xy = np.zeros([2,len(x)])
    xy[0] = x
    xy[1] = y
    xy = np.matrix(xy)

    MA = RotMatrix(theta)
    v = MA*xy
    
    v = np.array(v)
    
    xx = v[0]
    yy = v[1]
    
    return xx,yy
    
def disp_xy(x,y,disp_prec = 1e-3):
    
    tol = disp_prec/2
    
    mean_x = np.nanmean(x)
    mean_y = np.nanmean(y)
    std_x = np.nanstd(x)
    std_y = np.nanstd(y)
    xn = x - mean_x
    yn = y - mean_y
    
    # Defined in a 2-delta region
    N = np.int(np.ceil(std_x*4/disp_prec))
    M = np.int(np.ceil(std_y*4/disp_prec))
    C = CP.center(M,N)
    mg = MP.meshgridMN(M,N,xprec=disp_prec,yprec=disp_prec,C=C)
    
    Ixy = np.zeros([M,N])
    for ii in range(len(xn)):
        indxX = np.where(np.abs(mg.meshx-xn[ii])<tol)
        indxY = np.where(np.abs(mg.meshy-yn[ii])<tol)
        
        if len(indxX[0]) == 0 or len(indxY[0]) == 0: # Empty, out of boundaries
            continue
        else:
            ii = indxY[0][0]
            jj = indxX[1][0]
            Ixy[ii,jj] = Ixy[ii,jj] + 1
            
    return Ixy,mg