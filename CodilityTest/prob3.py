# Lattice next position problem
# Author: Yang Yu
# Date: 2018/10/16
# Description: Find the next position after due turn

import numpy as np

class line:
    a = []
    b = []
    c = []
    linedir = []
    
    def __init__(self,a=[],b=[],c=[]):
        self.a = a
        self.b = b
        self.c = c
    
    def lineconst(self,pt1,pt2):
        if pt1[0] == pt2[0]:
            self.a = 1
            self.b = 0
            self.c = -pt1[1]
        else:
            self.b = 1
            self.a = -(pt2[1]-pt1[1])/(pt2[0]-pt1[0])
            self.c = -self.a*pt1[0]-self.b*pt1[1]
            
        self.linedir = vec_normalize(pt2-pt1)
    
    def chkpt(self,pt,tol=1e-5):
        residual = np.abs(self.a*pt[0]+self.b*pt[1]+self.c)
        ptonlineflag = 0
        if residual < tol:
            ptonlineflag = 1
        
        return ptonlineflag
    
    def pendlineconst(self,pt):
        ptonlineflag = self.chkpt(pt)
        if ptonlineflag == 0:
            return -1
        else:
            if len(self.linedir) == 0: # Empty
                pendline = line()
                pendline.a = self.b
                pendline.b = -self.a
                pendline.c = -(pendline.a*pt[0]+pendline.b*pt[1])
            else:
                pendlinedir = np.zeros(2)
                pendlinedir[0] = self.linedir[1]
                pendlinedir[1] = -self.linedir[0]
                
                pendline = line()
                pendline.a = -pendlinedir[1]
                pendline.b = pendlinedir[0]
                pendline.c = -(pendline.a*pt[0]+pendline.b*pt[1])
                pendline.linedir = pendlinedir
            
            return pendline

def vec_normalize(v):
    v = v/np.sqrt(np.sum(v**2))
    
    return v

def integerchk(a):
    a_int = np.round(a)
    intflag = 0
    if np.abs(a-a_int)<1e-5:
        intflag = 1
    
    return intflag

def find_next_lattice_pt(Ax,Ay,Bx,By):
    
    Nextpt = np.zeros(2)
    
    pt1 = np.array([Ax,Ay])
    pt2 = np.array([Bx,By]) # Formatter
    
    dirAB = pt2-pt1
    
    lineAB = line()
    lineAB.lineconst(pt1,pt2)
    lineBC = lineAB.pendlineconst(pt2)
    
    # Find the next point according to the direction of lineBC
    BC_dir = lineBC.linedir
    # Calculate the x direction next point
    if np.abs(BC_dir[0] == 0 and BC_dir[1] == 0):
        return -1
    elif np.abs(BC_dir[0] == 0): # Search y direction
        increment = np.sign(BC_dir[1])
        pt_x = pt2[0]
        pt_y = pt2[1]+increment
    elif np.abs(BC_dir[1] == 0): # Search x direction
        increment = np.sign(BC_dir[0])
        pt_x = pt2[0]+increment
        pt_y = pt2[1]
    elif np.abs(BC_dir[0]) < np.abs(BC_dir[1]): # Search x direction
        increment = np.sign(BC_dir[0])
        pt_x = increment+pt2[0]
        while 1:
            pt_y = -(lineBC.c+lineBC.a*pt_x)/lineBC.b
            if integerchk(pt_y) == 1: # Not interger
                break
            else:
                pt_x = pt_x+increment
    elif np.abs(BC_dir[0]) >= np.abs(BC_dir[1]): # Search y direction
        increment = np.sign(BC_dir[1])
        pt_y = increment+pt2[1]
        while 1:
            pt_x = -(lineBC.c+lineBC.b*pt_y)/lineBC.a
            if integerchk(pt_x) == 1: # Not interger
                break
            else:
                pt_y = pt_y+increment
    
    Nextpt[0] = pt_x
    Nextpt[1] = pt_y
    
    return Nextpt
    
## Test routine
Ax = -1
Ay = 3
Bx = 3
By = 1

Nextpt = find_next_lattice_pt(Ax,Ay,Bx,By)
print(Nextpt)

Ax = 2
Ay = 2
Bx = 2
By = -3

Nextpt = find_next_lattice_pt(Ax,Ay,Bx,By)
print(Nextpt)