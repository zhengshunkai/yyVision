## Gettemplate from (ii,jj) position in I
# Author: Yang Yu
# Date: 2017/12/28
# Description: Get template from image

import numpy as np # For array operation
from meshop import meshgridi
from centerop import imcenter
    
def gettemplate_old(I,ii,jj,MT,NT): # Integer version
    
    # Get image size
    M,N = I.shape
    
    # Initialization of template: nans
    Itmp = np.empty((MT,NT))
    Itmp[:] = np.nan
    
    C = imcenter(Itmp)
    Itmp_mg = meshgridi(Itmp,C)
    
    # Calculate the position in the image I
    Itmp_mesh_mm = Itmp_mg.meshy + ii # y is vertical (ii)
    Itmp_mesh_nn = Itmp_mg.meshx + jj # x is horizontal (jj)
    
    for rr in range(MT):
        for ss in range(NT):
            mm = Itmp_mesh_mm[rr,ss]
            nn = Itmp_mesh_nn[rr,ss]
            
            if mm<0 or mm>=M or nn<0 or nn>=N: # Out of the boundaries
                continue
            else:
                Itmp[rr,ss] = I[mm,nn]
                
    return Itmp
    
def gettemplate(I,ii,jj,MT,NT): # Optimized version
    
    # Get image size
    M,N = I.shape
    
    # Initialization of template: nans
    Itmp = np.empty((MT,NT))
    Itmp[:] = np.nan
    
    C = imcenter(Itmp)
    Itmp_mg = meshgridi(Itmp,C)
    
    # Boundary check
    Itmp_mesh_mm = np.array(Itmp_mg.meshy[:,0] + ii,dtype='int') # y is vertical (ii)
    Itmp_mesh_nn = np.array(Itmp_mg.meshx[0,:] + jj,dtype='int') # x is horizontal (jj)
    
    ii_indx = np.where((Itmp_mesh_mm>=0)*(Itmp_mesh_mm<M))
    jj_indx = np.where((Itmp_mesh_nn>=0)*(Itmp_mesh_nn<N))
    
    Itmp[ii_indx[0][0]:ii_indx[0][-1]+1,jj_indx[0][0]:jj_indx[0][-1]+1] = I[Itmp_mesh_mm[ii_indx[0][0]]:Itmp_mesh_mm[ii_indx[0][-1]]+1,Itmp_mesh_nn[jj_indx[0][0]]:Itmp_mesh_nn[jj_indx[0][-1]]+1]
    
    return Itmp