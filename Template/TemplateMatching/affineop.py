
# Affine transformation
# Author: Yang Yu
# Date: 2018/4/10
# Description: NA

import numpy as np

def vecaffine(v,A=np.array([[1,0],[0,1]]),b=np.array([0,0])):
    v_e = np.matrix([v[0],v[1],1])
    A_e = np.matrix([[A[0,0],A[0,1],b[0]],[A[1,0],A[1,1],b[1]],[0,0,1]])
    y_e = A_e*v_e.T
    
    y_e_a = np.array(y_e.T)
    y = y_e_a[0][0:2]
    
    return y

def affineRotA(theta,radflag = 1): # rad
    if radflag == 0:
        theta = theta/180*pi
    
    RotA = np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    
    return RotA
    
def vecaffine_mult(v_mult,A=np.array([[1,0],[0,1]]),b=np.array([0,0])): # v_mult is M lines 2 cols
    _,N = v_mult.shape
    y_mult = np.zeros(v_mult.shape)
    for ii in range(N):
        y_mult[:,ii] = vecaffine(v_mult[:,ii],A,b)
    
    return y_mult
    
def CtrlPt_pos_deneg(pv):
    x = pv[1,:]
    y = pv[0,:]
    x_min = np.nanmin(x)
    y_min = np.nanmin(y)
    x_r = x - x_min
    y_r = y - y_min
    
    pv_r = np.zeros(pv.shape)
    pv_r[0,:] = y_r
    pv_r[1,:] = x_r
    
    return pv_r
    
# Generate M unique number from a range
def randindxgen(M,str=0,fin=255):
    
    candidates = np.array(range(str,fin))
    # Shuffle the candidates
    np.random.shuffle(candidates)
    UniqRandomIndx = candidates[0:M]
    
    return UniqRandomIndx
    
def shift_loc_pts(CtrlPt_local_pos,T=0.05):
    delta = T/2/np.sqrt(2) - 1e-7
    M,N = CtrlPt_local_pos.shape
    
    shift_array = np.zeros(CtrlPt_local_pos.shape)
    shift_array_x = np.random.rand(M,N)[0]*delta
    shift_array_y = np.random.rand(M,N)[0]*delta
    shift_array[1,:] = shift_array_x
    shift_array[0,:] = shift_array_y
    CtrlPt_local_pos = CtrlPt_local_pos + shift_array
        
    return CtrlPt_local_pos