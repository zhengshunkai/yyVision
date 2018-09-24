## Differentiate routine
# Author: Yang Yu
# Date: 2018/4/27
# Description: dy/dx

import math
from math import pi
import numpy as np

import matplotlib.pyplot as plt

def deline2(x,y):
    
    deri = np.zeros(len(x))
    for ii in range(len(x)):
        if ii == 0:
            delta_x = x[ii+1]-x[ii]
            delta_y = y[ii+1]-y[ii]
        elif ii == len(x)-1:
            delta_x = x[ii]-x[ii-1]
            delta_y = y[ii]-y[ii-1]
        else:
            delta_x = x[ii+1] - x[ii-1]
            delta_y = y[ii+1] - y[ii-1]
    
        deri[ii] = delta_y/delta_x
    
    return deri
    
def deline_l1(x,y):
    
    deri = np.zeros(len(x))
    deri[0] = np.nan
    for ii in range(1,len(x)):
        delta_x = x[ii]-x[ii-1]
        delta_y = y[ii]-y[ii-1]
            
        deri[ii] = delta_y/delta_x
    
    return deri

def deline_r1(x,y):
    
    deri = np.zeros(len(x))
    deri[-1] = np.nan
    for ii in range(len(x)-1):
        delta_x = x[ii]-x[ii-1]
        delta_y = y[ii]-y[ii-1]
    
        deri[ii] = delta_y/delta_x
        
    return deri
    
def diff1D(x):
    
    x_1= np.zeros(len(x))
    x_1[0:len(x)-1] = x[1:len(x)]
    x_1[-1] = x[-2]
    
    dx = x_1-x
    dx[-1] = dx[-2]
    
    return dx