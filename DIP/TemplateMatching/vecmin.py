## Vec/array operations
# Author: Yang Yu
# Date: 2018/4/10
# Description: NA

import numpy as np

def nanargmaxm(A,N=1):
    tmp = np.float32(A.copy())
    tmp = np.reshape(tmp,[1,-1])
    tmp = tmp[0]
    for ii in range(N):
        indx = np.nanargmax(tmp)
        tmp[indx] = -np.inf
    tmp2 = np.reshape(tmp,A.shape)
    indx = np.where(tmp2 == -np.inf)
    
    return indx
    
def nanmaxm(A,N=1):
    tmp = np.float32(A.copy())
    tmp = np.reshape(tmp,[1,-1])
    tmp = tmp[0]
    for ii in range(N):
        indx = np.nanargmax(tmp)
        tmp[indx] = -np.inf
    tmp2 = np.reshape(tmp,A.shape)
    indx = np.where(tmp2 == -np.inf)
    
    val = A[indx[0],indx[1]]
    
    return val
    
def vecmin(v):
    
    vmin = np.array(list(set(v.tolist())))
    
    return vmin
    
def vecmin2(vec_ens):
    
    M,N = vec_ens.shape
    
    vec_ens_l = vec_ens.tolist()
    vec_ens_min_l = []
    while 1:
        vectmp = vec_ens_l[0] # Get the first element
        
        indx = findvec(vec_ens_l,vectmp)
        
        # Search will find at least 1 element @ pos 0
        vec_ens_min_l.append(vectmp)
        for jj in range(len(indx[0])):
            vec_ens_l.remove(vectmp) # Remove all vectmp
        
        if len(vec_ens_l) == 0: # empty
            break
    
    vec_ens_min = np.array(vec_ens_min_l)
    
    return vec_ens_min
    
def findvec(vec_ens,vec,tol=1e-3):
    
    vec_ens = np.array(vec_ens)
    vec = np.array(vec) # In case list is used
    indx = np.where((np.abs(vec_ens[:,0]-vec[0]) < tol)*(np.abs(vec_ens[:,1]-vec[1]) < tol))
    
    return indx