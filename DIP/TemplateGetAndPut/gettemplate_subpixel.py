## Gettemplate from (ii,jj) position in I
# Author: Yang Yu
# Date: 2017/12/28
# Description: Get template from image

import numpy as np
import scipy as sp
from meshop import meshgridi
from centerop import imcenter
                   
def gettemplatesb(I,ii,jj,MT,NT,method='cubic'): # Subpixel version

    # Get image size
    M,N = I.shape
    
    # Initialization of template: nans
    Itmp = np.empty((MT,NT))
    Itmp[:] = np.nan

    C = imcenter(Itmp)
    Itmp_mg = meshgridi(Itmp,C)    
    I_mg = meshgridi(I)
    I_mesh_x = I_mg.meshx[1,:]
    I_mesh_y = I_mg.meshy[:,1]
    
    # Calculate the position in the image I
    # ii, jj subpixel
    Itmp_mesh_mm = Itmp_mg.meshy + ii # y is vertical (ii)
    Itmp_mesh_nn = Itmp_mg.meshx + jj # x is horizontal (jj)
    Itmp_mesh_xx = Itmp_mesh_nn[1,:]
    Itmp_mesh_yy = Itmp_mesh_mm[:,1]
    
    # Interpolation
    f_fit = sp.interpolate.interp2d(I_mesh_x,I_mesh_y,I,kind=method)
    Itmp = f_fit(Itmp_mesh_xx,Itmp_mesh_yy)
                
    return Itmp

