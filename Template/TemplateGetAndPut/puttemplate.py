import numpy as np # For array operation
from meshop import meshgridi
from centerop import imcenter

def puttemplate(I,tmp,ptx,pty):
    
    M,N = I.shape
    
    # Boundary check
    C = imcenter(tmp)
    tmp_mg = meshgridi(tmp,C)
    
    # Boundary check
    Itmp_mesh_mm = np.array(tmp_mg.meshy[:,0] + pty,dtype='int') # y is vertical (ii)
    Itmp_mesh_nn = np.array(tmp_mg.meshx[0,:] + ptx,dtype='int') # x is vertical (jj)
    
    ii_indx = np.where((Itmp_mesh_mm>=0)*(Itmp_mesh_mm<M))
    jj_indx = np.where((Itmp_mesh_nn>=0)*(Itmp_mesh_nn<N))
    
    I[Itmp_mesh_mm[ii_indx[0][0]]:Itmp_mesh_mm[ii_indx[0][-1]]+1,Itmp_mesh_nn[jj_indx[0][0]]:Itmp_mesh_nn[jj_indx[0][-1]]+1] = tmp[ii_indx[0][0]:ii_indx[0][-1]+1,jj_indx[0][0]:jj_indx[0][-1]+1]
    
    return I