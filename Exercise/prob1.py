# Stable period count
# Author: Yang Yu
# Description: Count the stable period within a given discrete signal
# Use numpy (numpy)

import numpy as np

def count_stable(A,Amin=-1e9,Amax=1e9): 
    
    p_num = 0 # Number of stable period
    p = [] # Stable periods record
    
    L = len(A) # Find the length of A
    p_increment = -1 # Period increment
    p_str = -1 # Start of current period
    p_end = -1 # End of current period
    
    # Value check
    A = np.array(A,dtype='float')
    indx = np.where((A<Amin)+(A>Amax))
    A[indx] = np.nan
    
    for ii in range(0,L):
        
        p_str = ii
        
        if p_str == L-1: # Last element received
            break
        else:
            for jj in range(p_str+1,L):
                if jj == p_str+1: # The next element
                    p_end = jj
                    p_increment = A[p_str+1]-A[p_str]
                else:
                    increment_tmp = A[jj]-A[jj-1]
                    if increment_tmp == p_increment:
                        p_end = jj
                        if jj >= p_str+2: # minimum length 2
                            p.append([p_str,p_end])
                            p_num = p_num + 1
                            
                            if p_num > 1e9:
                                return -1
                    else:
                            break
    return p_num,p
                    
## Test routine
A = np.array([-1,1,3,3,3,2,3,2,1,0])
p_num,p = count_stable(A)
