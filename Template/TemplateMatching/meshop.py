## Generate a meshgrid for an image
# Author: Yang Yu
# Date: 2017/12/28
# Description: meshgrid generation

def meshgrid(I,C=[0,0],sparseflag=0):
    import numpy as np
    import math
    # Class definition
    class mg:
        meshx = []
        meshy = []
        meshr = []
        xaxis = []
        yaxis = []
        xprec = []
        yprec = []
    
    Vord,Hord = I.shape
    xprec = 1/(Hord-1)
    yprec = 1/(Vord-1)
    xaxis = np.arange(0,Hord,1)*xprec # x axis is horiztonal
    yaxis = np.arange(0,Vord,1)*yprec # y axis is Vertical
    xaxis = xaxis - C[1]*xprec
    yaxis = yaxis - C[0]*yprec
    
    if sparseflag == 0:
        meshx,meshy = np.meshgrid(xaxis,yaxis)
    elif sparseflag == 1:
        meshx,meshy = np.meshgrid(xaxis,yaxis,sparse=True) # Using sparse matrix
    meshr = (meshx**2 + meshy**2)**0.5
    
    mg.xprec = xprec
    mg.yprec = yprec
    mg.xaxis = xaxis
    mg.yaxis = yaxis
    mg.meshx = meshx
    mg.meshy = meshy
    mg.meshr = meshr
    
    return mg
    
def meshgridi(I,C=[0,0],sparseflag=0):
    import numpy as np
    import math
    # Class definition
    class mg:
        meshx = []
        meshy = []
        meshr = []
        xaxis = []
        yaxis = []
        xprec = []
        yprec = []
    
    Vord,Hord = I.shape
    xaxis = np.arange(0,Hord,1) # x axis is horiztonal
    yaxis = np.arange(0,Vord,1) # y axis is Vertical
    xaxis = xaxis - C[1]
    yaxis = yaxis - C[0]
    
    if sparseflag == 0:
        meshx,meshy = np.meshgrid(xaxis,yaxis)
    elif sparseflag == 1:
        meshx,meshy = np.meshgrid(xaxis,yaxis,sparse=True) # Using sparse matrix
    meshr = (meshx**2 + meshy**2)**0.5
    
    mg.xprec = 1
    mg.yprec = 1
    mg.xaxis = xaxis
    mg.yaxis = yaxis
    mg.meshx = meshx
    mg.meshy = meshy
    mg.meshr = meshr
    
    return mg
    
def meshgridMN(Vord,Hord,xprec=1,yprec=1,C=[0,0],Shift=[0,0],sparseflag=0):
    import numpy as np
    import math
    # Class definition
    class mg:
        meshx = []
        meshy = []
        meshr = []
        xaxis = []
        yaxis = []
        xprec = []
        yprec = []
    
    xaxis = np.arange(0,Hord,1)*xprec # x axis is horiztonal
    yaxis = np.arange(0,Vord,1)*yprec # y axis is Vertical
    xaxis = xaxis-C[1]*xprec+Shift[1]
    yaxis = yaxis-C[0]*yprec+Shift[0]
    
    if sparseflag == 0:
        meshx,meshy = np.meshgrid(xaxis,yaxis)
    elif sparseflag == 1:
        meshx,meshy = np.meshgrid(xaxis,yaxis,sparse=True) # Using sparse matrix
    meshr = (meshx**2 + meshy**2)**0.5
    
    mg.xprec = xprec
    mg.yprec = yprec
    mg.xaxis = xaxis
    mg.yaxis = yaxis
    mg.meshx = meshx
    mg.meshy = meshy
    
    return mg