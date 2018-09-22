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