## Gettemplate optimization version
# Author: Yang Yu
# Date: 2018/4/28
# Description: NA

import math
from math import pi
import numpy as np
import time

import gettemplate as GT
import puttemplate as PT
import scipy.ndimage.interpolation as interp

import matplotlib.pyplot as plt

I = np.reshape(np.arange(100000),[1000,100])

MT = 100
NT = 100
ii = 1
jj = 1

# Routine 1
tic = time.clock()
tmp1 = GT.gettemplate_old(I,ii,jj,MT,NT)
toc = time.clock()
elapsetime1 = toc - tic

# Routine 2
tic = time.clock()
tmp2 = GT.gettemplate(I,ii,jj,MT,NT)
toc = time.clock()
elapsetime2 = toc - tic

print(elapsetime1,elapsetime2,elapsetime1/elapsetime2)

plt.figure()
plt.imshow(I)
plt.figure()
plt.subplot(1,2,1)
plt.imshow(tmp1)
plt.subplot(1,2,2)
plt.imshow(tmp2)

# Puttemplate test
tmp = np.random.randint(0,255,[100,100])
tmp_rot = interp.rotate(tmp,45/pi*180,reshape=1)

plt.figure()
plt.imshow(tmp)
plt.figure()
plt.imshow(tmp_rot)

I = np.zeros([1000,1000])
pt = [50,50]
J = PT.puttemplate(I,tmp,pt[1],pt[0])

plt.figure()
plt.imshow(J)

plt.show()
