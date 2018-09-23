## Combination operations
# Author: Yang Yu
# Date: 2018/4/13
# Description: https://stackoverflow.com/questions/16003217/n-d-version-of-itertools-combinations-in-numpy

import numpy as np
from itertools import combinations, chain
from scipy.special import comb

def comb_index(n, k):
    count = comb(n, k, exact=True)
    index = np.fromiter(chain.from_iterable(combinations(range(n), k)), 
                        int, count=count*k)
    return index.reshape(-1, k)
    
def comb_num(n, k):
    
    upper_part = 1
    down_part = 1
    for ii in range(k):
        upper_part = upper_part * (n-ii)
        down_part = down_part * (ii+1)
    
    combnum = np.int(upper_part/down_part)
    
    return combnum
    
def combinations(v,k):
    v_num = len(v)
    
    comb = []
    combnum = 0
    
    if v_num >= k:
        combnum = CI.comb_num(v_num,k)
        comb = np.zeros([combnum,k]) # initialization
        comb_index = CI.comb_index(v_num,k)
        M,N = comb_index.shape
        
        for ii in range(M):
            for jj in range(k):
                comb[ii,jj] = v[comb_index[ii,jj]]
        
    return comb,combnum