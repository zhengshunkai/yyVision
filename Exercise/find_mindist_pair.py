# Find minimum distance pair
# Author: Yang Yu
# Date: 2018/10/16
# Description: Find minimum distance pair within an array
# Find all minimum distance pair

import numpy as np

def find_mindist_pair(A):
    
    minval = -np.inf # minimum value
    minpos = -1 # minimum pair position
    
    L = len(A)
    A_shifted = np.zeros(A.shape)
    A_shifted[0:L-1] = A[1:L]
    
    A_diff = np.abs(A_shifted-A)
    A_diff[-1] = np.inf # Set the final element to infinite
    
    minval = np.nanmin(A_diff)
    minpos = np.where(A_diff == minval)
    
    return minval,minpos


## Test routine
A = np.array([0,3,3,7,5,3,11,1],dtype='int32')
# -2,147,483,648 32 bits
minval,minpos = find_mindist_pair(A)

print(minval)
print(minpos)