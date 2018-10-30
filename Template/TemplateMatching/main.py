## Template matching: a control points & ensemble algorithm approach
# Author: Yang Yu
# Date: 2018/4/12
# Description: Outliner handling via an ensemble algorithm approach

import math
from math import pi
import numpy as np
import scipy as sp
import scipy.signal as spsg
import cv2
import time

import centerop
import meshop
import disttable as DT
import affineop as AF
import vecmin
from TemplateMatching import control_points_template_match_outH,control_points_template_match
from comb_index import comb_index

import matplotlib.pyplot as plt

## Functions
# 

## Global map generation
# Global map/dist table generation
GMconfig = DT.globalMap_config()
GMconfig.CtrlPtNum = 100

# Local map configuration
local_length = 30 # meter
local_width = 30 # meter
local_center = centerop.center(local_width,local_length)

global_length = GMconfig.global_length
global_width = GMconfig.global_width
simu_prec = GMconfig.simu_prec
disp_prec = GMconfig.disp_prec
CtrlPtNum = GMconfig.CtrlPtNum

newGlobalMapflag = 0
xmlfile = 'GlobalMap.xml'
xmlfile_Ext = 'GlobalMap_Ext.xml'

tol = 0.1*2*np.sqrt(2)
angtol = 5/180*pi

if newGlobalMapflag == 1:
    # Generate the CtrlPt positions (x,y unit meter)
    CtrlPt_x = np.random.randint(0,np.uint(global_length/simu_prec)-1,size=CtrlPtNum)*simu_prec
    CtrlPt_y = np.random.randint(0,np.uint(global_width/simu_prec)-1,size=CtrlPtNum)*simu_prec
    CtrlPt_pos = np.zeros([2,CtrlPtNum],dtype='float')
    CtrlPt_pos[1,:] = CtrlPt_x # x->cols
    CtrlPt_pos[0,:] = CtrlPt_y # y->Lines
    
    # Display CtrlPt position/GlobalMap demonstration
    CtrlPt_pos_disp,GlobalMap = DT.dispGlobalMap(CtrlPt_pos,GMconfig)
    # GlobalMap2 = DT.dispMap(CtrlPt_pos,GMconfig)
    
    # Sorting CtrlPts
    CtrlPt_pos_unsorted = CtrlPt_pos.copy()
    CtrlPt_pos_sorted = DT.indx_sort(CtrlPt_pos_unsorted,GMconfig)
    
    # Record distances between all points
    global_dist_table = DT.dist_table_gen(CtrlPt_pos_sorted)
    global_ang_table = DT.ang_table_gen(CtrlPt_pos_sorted)
    
    # Temporary save
    # Creat filestorage, store only the dist table + CtrlPt position list
    fs = cv2.FileStorage(xmlfile,cv2.FILE_STORAGE_WRITE)
    fs.write('ref_pos_sorted',CtrlPt_pos_sorted)
    fs.write('global_dist_table',global_dist_table)
    fs.write('global_ang_table',global_ang_table)
    fs.release()
    
    # Creat filestorage
    fs = cv2.FileStorage(xmlfile_Ext,cv2.FILE_STORAGE_WRITE)
    fs.write('Ref_pos',CtrlPt_pos)
    fs.write('Ref_pos_disp',CtrlPt_pos_disp)
    fs.write('GlobalMap',GlobalMap)
    fs.write('ref_pos_sorted',CtrlPt_pos_sorted)
    fs.write('global_dist_table',global_dist_table)
    fs.write('global_ang_table',global_ang_table)
    fs.release()
else:
    # Read the filestorage
    fs = cv2.FileStorage(xmlfile,cv2.FILE_STORAGE_READ)
    fn = fs.getNode('ref_pos_sorted')
    CtrlPt_pos_sorted = fn.mat()
    fn = fs.getNode('global_dist_table')
    global_dist_table = fn.mat()
    fn = fs.getNode('global_ang_table')
    global_ang_table = fn.mat()
    fs.release()
    
    fs = cv2.FileStorage(xmlfile_Ext,cv2.FILE_STORAGE_READ)
    fn = fs.getNode('GlobalMap')
    GlobalMap = fn.mat()
        
plt.figure()
plt.imshow(GlobalMap)
plt.figure()
plt.imshow(global_dist_table)
plt.figure()
plt.imshow(global_ang_table)

## Control points template matching simu
# A sweeping approach
sweep_step_x = 10 # meter
sweep_step_y = 10 # meter
x_margin = 5 # meter
y_margin = 5 # meter margin areas not swiped

sweep_x_num = np.floor((global_length-2*x_margin)/sweep_step_x)
sweep_y_num = np.floor((global_width-2*y_margin)/sweep_step_y)
sweep_mg = meshop.meshgridMN(sweep_y_num,sweep_x_num,sweep_step_y,sweep_step_x,Shift=[y_margin,x_margin])

selindx = 4
for ii in range(selindx-1,selindx):
    for jj in range(selindx-1,selindx):
        
        ## Generate the local template
        # Get the template from the CtrlPt_pos
        center_x = sweep_mg.meshx[ii,jj]
        center_y = sweep_mg.meshy[ii,jj]
        win_x_half = local_length/2
        win_y_half = local_width/2
        
        # Get the local map
        CtrlPt_x = CtrlPt_pos_sorted[1,:]
        CtrlPt_y = CtrlPt_pos_sorted[0,:]
        indx = np.where((CtrlPt_x <= (center_x+win_x_half))*(CtrlPt_x >= (center_x-win_x_half))*(CtrlPt_y <= (center_y+win_y_half))*(CtrlPt_y >= (center_y-win_y_half)))
        CtrlPt_local_pos = CtrlPt_pos_sorted[:,indx[0]]
        localmap = DT.dispMap(CtrlPt_local_pos,GMconfig)
        local_dist_table = DT.dist_table_gen(CtrlPt_local_pos)
        local_ang_table = DT.ang_table_gen(CtrlPt_pos_sorted)
        
        plt.figure()
        plt.imshow(localmap)
        plt.figure()
        plt.imshow(local_dist_table)
        plt.figure()
        plt.imshow(local_ang_table)
        
        ## Local CtrlPt pos distortion
        CtrlPt_local_pos_distorted = CtrlPt_local_pos.copy() # No distortion case
        
        # Occlusion
        M = 7
        local_CtrlPt_num_no_occlusion = len(CtrlPt_local_pos_distorted[0,:])
        randindx = AF.randindxgen(M,0,local_CtrlPt_num_no_occlusion)
        # randindx = np.array([17, 21,  2,  0])
        CtrlPt_local_pos_distorted = CtrlPt_local_pos_distorted[:,randindx]
        # CtrlPt_local_pos_distorted = CtrlPt_local_pos_distorted.copy()
        
        # Introduce noise
        CtrlPt_local_pos_distorted = AF.shift_loc_pts(CtrlPt_local_pos_distorted,T=tol)
        
        # Add false CtrlPt
        N = 2
        false_CtrlPt_x = np.random.randint(0,np.uint(local_length/simu_prec)-1,size=N)*simu_prec + (center_x-win_x_half)
        false_CtrlPt_y = np.random.randint(0,np.uint(local_width/simu_prec)-1,size=N)*simu_prec + (center_y-win_y_half)
        tmplist = CtrlPt_local_pos_distorted[1,:].tolist()
        tmplist = tmplist + false_CtrlPt_x.tolist()
        CtrlPt_local_pos_distorted_x = np.array(tmplist)
        tmplist = CtrlPt_local_pos_distorted[0,:].tolist()
        tmplist = tmplist + false_CtrlPt_y.tolist()
        CtrlPt_local_pos_distorted_y = np.array(tmplist)
        CtrlPt_local_pos_distorted = np.zeros([2,M+N],dtype='float')
        CtrlPt_local_pos_distorted[1,:] = CtrlPt_local_pos_distorted_x
        CtrlPt_local_pos_distorted[0,:] = CtrlPt_local_pos_distorted_y
        
        # Rotation
        theta = 45/180*pi
        A = AF.affineRotA(theta)
        CtrlPt_local_pos_distorted = AF.vecaffine_mult(CtrlPt_local_pos_distorted,A)
        
        # Negative value handling
        CtrlPt_local_pos_distorted = AF.CtrlPt_pos_deneg(CtrlPt_local_pos_distorted)
        
        localmap_distorted = DT.dispMap(CtrlPt_local_pos_distorted,GMconfig)
        local_dist_table_distorted = DT.dist_table_gen(CtrlPt_local_pos_distorted)
        local_ang_table_distorted = DT.ang_table_gen(CtrlPt_local_pos_distorted)
        
        plt.figure()
        plt.imshow(localmap_distorted)
        plt.figure()
        plt.imshow(local_dist_table_distorted)
        
        ## Control points template matching
        start_time = time.clock()
        # hit,local_CtrlPt_sel,search_candidates,search_candidates_prev = control_points_template_match_outH(local_dist_table_distorted,global_dist_table,T=tol) # Find the control points in the global map (outliner handling version)
        hit,local_CtrlPt_sel,search_candidates = control_points_template_match(local_dist_table_distorted,local_ang_table_distorted,CtrlPt_local_pos_distorted,global_dist_table,global_ang_table,disttol=tol,angtol=angtol)
        stop_time = time.clock()
        elapsetime = stop_time-start_time
        print(elapsetime)
        
        if not hit: 
            if not search_candidates:
                print("No findings.")
            else:
                print("No unique hit")
                SC_localMap = DT.dispMap(CtrlPt_local_pos_distorted[:,local_CtrlPt_sel], GMconfig, step=1)
                plt.figure()
                plt.imshow(SC_localMap)

                for candindx in range(len(search_candidates)):
                    SC_CtrlPt_pos = CtrlPt_pos_sorted[:,search_candidates[candindx].ctrl_pts]
                    SC_CtrlPt_pos_disp,SC_GlobalMap = DT.dispGlobalMap(SC_CtrlPt_pos,GMconfig,step=1)
                    
                    plt.figure()
                    plt.imshow(SC_GlobalMap)
        else:
            print("Hit")
            hit_CtrlPt_pos = CtrlPt_pos_sorted[:,hit.ctrl_pts]
            hit_CtrlPt_pos_disp,hit_GlobalMap = DT.dispGlobalMap(hit_CtrlPt_pos,GMconfig,step=1)
            hit_localMap = DT.dispMap(CtrlPt_local_pos_distorted[:,local_CtrlPt_sel], GMconfig, step=1)
            
            plt.figure()
            plt.imshow(hit_localMap)
            plt.figure()
            plt.imshow(hit_GlobalMap)
                
## Debug show plots
plt.show()