## Find the center of an image
# Author: Yang Yu
# Date: 2017/12/28
# Description: find the center of an image, image I is a matrix or 2D array
# Return a 2 by 1 matrix

def imcenter(I,origflag=0):
    import math
    import numpy as np
    import scipy as sp
    
    M,N = I.shape
    C = center(M,N,origflag)

    
    return C
    
def center(M,N,origflag=0): 
    import math
    C = [math.floor(M/2)+1,math.floor(N/2)+1]
    
    if origflag == 0:
        C = [C[0]-1,C[1]-1] 
    
    return C