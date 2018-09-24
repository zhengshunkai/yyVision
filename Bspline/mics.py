## Differentiate routine
# Author: Yang Yu
# Date: 2018/4/27
# Description: dy/dx

import math
from math import pi
import numpy as np

import matplotlib.pyplot as plt

def deline(line): # dy/dx
    
    x = np.array(line)[:,0]
    y = np.array(line)[:,1]
    
    deri = np.zeros(line.shape)
    deri[:,0] = x.copy()
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
        
        deri[ii,1] = delta_y/delta_x
    
    return deri
    