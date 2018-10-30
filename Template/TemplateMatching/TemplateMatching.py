## Control point template matching special case
# Author: Yang Yu
# Date: 2018/4/8
# Description: NA

import math
from math import pi
import numpy as np
import scipy as sp
import scipy.signal as spsg
import cv2

import vecmin
import disttable as DT
from comb_index import comb_index

import matplotlib.pyplot as plt

class searchcand: # Initially want to cover all Topo structure description (generic)
    ctrl_pts = []
    ctrl_pts_minvec = []
    local_pts = []
    
    def __init__(self,ctrl_pts=np.array([])):
        self.ctrl_pts = ctrl_pts
        if not ctrl_pts.tolist():
            self.ctrl_pts_minvec = []
        else:
            self.set_ctrl_pts(ctrl_pts=ctrl_pts)
            
    def set_local_pts(self,local_pts):
        self.local_pts = local_pts
    
    def set_ctrl_pts(self,ctrl_pts):
        self.ctrl_pts_minvec = vecmin.vecmin(ctrl_pts)

## Functions
def dist_table_shriek(dist_table):
    # Keep only upper half
    dist_table_half = np.zeros(dist_table.shape)
    M,N = dist_table.shape
    
    for ii in range(M):
        for jj in range(N):
            if ii<jj:
                dist_table_half[ii,jj] = dist_table[ii,jj]
    
    return dist_table_half

def ctrl_pts2dist_table(ctrl_pts,dist_table):
    
    ctrl_pts_dist_table_lines = dist_table[ctrl_pts,:]
    ctrl_pts_dist_table = ctrl_pts_dist_table_lines[:,ctrl_pts]
            
    return ctrl_pts_dist_table

def add_local_ctrl_points(local_ctrl_pts,local_CtrlPtNum):
    
    if len(local_ctrl_pts) == local_CtrlPtNum:
        local_ctrl_pts_new = local_ctrl_pts
        return local_ctrl_pts_new
    else:
        local_ctrl_pts_new = np.zeros(len(local_ctrl_pts)+1,dtype='uint')
        local_ctrl_pts_new[0:len(local_ctrl_pts)] = local_ctrl_pts
    
    for ii in range(local_CtrlPtNum):
        indx = np.where(local_ctrl_pts == ii)
        if not indx[0].tolist(): # Empty
            local_ctrl_pts_new[-1] = ii
            break
        else:
            continue
    
    return local_ctrl_pts_new

def swap(a):
    tmp = a[0]
    a[0] = a[1]
    a[1] = tmp
    
    return a

def addNewCtrlPt(a,newpt):
    al = a.tolist()
    al.append(newpt)
    b = np.array(al)
    
    return b
    
## Check dist tables
def dist_table_check(local_ctrl_pts_dist_table,global_dist_table,T):
    
    local_ctrl_pts_num = len(local_ctrl_pts_dist_table[0,:])
    global_dist_table_half = dist_table_shriek(global_dist_table)
    
    for jj in range(1,local_ctrl_pts_num): # Excluding the [0,0] points, starts from [0,1]
        # Search routine
        # Use the first line for search and the columns to check (dist table)
        dist_target = local_ctrl_pts_dist_table[0,jj]
        
        # The first point, global search the segment that corresponds to max dist
        if jj == 1:
            
            search_candidates = []
            search_candidates_prev = []
            
            candpos = np.where(np.abs(global_dist_table_half-dist_target)<T)
            candnum = len(candpos[0]) # Determine candidates number
            
            for candindx in range(candnum):
                searchcand_ctrl_pts = np.array([candpos[0][candindx],candpos[1][candindx]])
                searchcand_tmp = searchcand(searchcand_ctrl_pts)
                search_candidates.append(searchcand_tmp)
                    
        # The second point, distinguish the order of the found control points
        # (P1,P2)->(p1,p2) or (P1,P2)->(p2,p1)
        elif jj == 2: 
            
            search_candidates_prev = search_candidates.copy() # previous level search candidates (1st points)
            search_candidates = []
            
            candnum = len(search_candidates_prev)
            for candindx in range(candnum): # A candidate by candidate search routine
                
                cand_ctrl_pts = search_candidates_prev[candindx].ctrl_pts
                
                # ord 0: (P1,P2)->(p1,p2)
                cand_search_ROI = global_dist_table[cand_ctrl_pts[0],:]
                new_ctrl_pts_pos_ensemble = np.where(np.abs(cand_search_ROI-dist_target)<T)
                new_ctrl_pts_ensemble = new_ctrl_pts_pos_ensemble[0] # 1D search
                
                # Check if the third ctrl point is valid
                if not new_ctrl_pts_ensemble.tolist():
                    pass
                else:
                    local_check_val = local_ctrl_pts_dist_table[1,2]
                    for kk in range(len(new_ctrl_pts_ensemble)):
                        global_check_val = global_dist_table[cand_ctrl_pts[1],new_ctrl_pts_ensemble[kk]]
                        if np.abs(global_check_val-local_check_val)<T:
                            searchcand_ctrl_pts = addNewCtrlPt(cand_ctrl_pts,new_ctrl_pts_ensemble[kk])
                            searchcand_tmp = searchcand(searchcand_ctrl_pts)
                            search_candidates.append(searchcand_tmp)
                    
                # ord 1: (P1,P2)->(p2,p1)
                swap(cand_ctrl_pts)
                cand_search_ROI = global_dist_table[cand_ctrl_pts[0],:]
                new_ctrl_pts_pos_ensemble = np.where(np.abs(cand_search_ROI-dist_target)<T)
                new_ctrl_pts_ensemble = new_ctrl_pts_pos_ensemble[0] # 1D search
                
                # Check if the third ctrl point is valid
                if not new_ctrl_pts_ensemble.tolist():
                    pass
                else:
                    local_check_val = local_ctrl_pts_dist_table[1,2]
                    for kk in range(len(new_ctrl_pts_ensemble)):
                        global_check_val = global_dist_table[cand_ctrl_pts[1],new_ctrl_pts_ensemble[kk]]
                        if np.abs(global_check_val-local_check_val)<T:
                            searchcand_ctrl_pts = addNewCtrlPt(cand_ctrl_pts,new_ctrl_pts_ensemble[kk])
                            searchcand_tmp = searchcand(searchcand_ctrl_pts)
                            search_candidates.append(searchcand_tmp)
                        
        # In this case, multiple checks are conducted           
        elif jj>2:
            
            search_candidates_prev = search_candidates.copy() # previous level search candidates (1st points)
            search_candidates = []
            
            candnum = len(search_candidates_prev)
            for candindx in range(candnum): # A candidate by candidate search routine
                
                cand_ctrl_pts = search_candidates_prev[candindx].ctrl_pts
                
                cand_search_ROI = global_dist_table[cand_ctrl_pts[0],:]
                new_ctrl_pts_pos_ensemble = np.where(np.abs(cand_search_ROI-dist_target)<T)
                new_ctrl_pts_ensemble = new_ctrl_pts_pos_ensemble[0] # 1D search
                
                # Check if the new ctrl points are valid
                if not new_ctrl_pts_ensemble.tolist():
                    pass
                else:
                    for kk in range(len(new_ctrl_pts_ensemble)):
                        checkflag = 1
                        for ll in range(1,jj):
                            local_check_val = local_ctrl_pts_dist_table[ll,jj]
                            global_check_val = global_dist_table[cand_ctrl_pts[ll],new_ctrl_pts_ensemble[kk]]
                            if np.abs(global_check_val-local_check_val)<T:
                                pass
                            else:
                                checkflag = 0
                                break
                    
                        if checkflag == 1:
                            searchcand_ctrl_pts = addNewCtrlPt(cand_ctrl_pts,new_ctrl_pts_ensemble[kk])
                            searchcand_tmp = searchcand(searchcand_ctrl_pts)
                            search_candidates.append(searchcand_tmp)
        
        candnum = len(search_candidates)
        if candnum == 0: # No candidates found
            break
        
    return search_candidates,search_candidates_prev

def ang_table_check(search_candidates,global_ang_table,local_ang_table,angtol = 5/180*pi):
    
    search_candidates_angchk = []
    for candindx in range(len(search_candidates)):
        ctrl_pts = search_candidates[candindx].ctrl_pts
        local_CtrlPt_sel = search_candidates[candindx].local_pts
        local_ang_table_sel = ctrl_pts2dist_table(local_CtrlPt_sel,local_ang_table)
        
        local_ang_table_norm = local_ang_table_sel - local_ang_table_sel[0,1]
        global_ang_table_norm = global_ang_table - global_ang_table[ctrl_pts[0],ctrl_pts[1]]
        
        M = len(ctrl_pts)
        chkflag = 1
        for ii in range(M):
            for jj in range(ii+1,M):
                localang = DT.angadjrange(local_ang_table_norm[ii,jj])
                globalang = DT.angadjrange(global_ang_table_norm[ctrl_pts[ii],ctrl_pts[jj]])
                if np.abs(globalang-localang) < angtol:
                    continue
                else:
                    chkflag = 0
                    break
            if chkflag == 0:
                break
        
        if chkflag == 1:
            search_candidates_angchk.append(search_candidates[candindx])
    
    return search_candidates_angchk

def control_points_template_match(local_dist_table,global_dist_table,T=0.05):
    
    # Ctrl_pts
    local_ctrl_pts = []
    
    ## initialization: find the initial local ctrl points topology (dist_table)
    # Find the maximum distance line/Pts pair
    local_CtrlPtNum,_ = local_dist_table.shape
    local_dist_table_half = dist_table_shriek(local_dist_table)
    
    maxdist = np.nanmax(local_dist_table)
    indx = np.where(local_dist_table_half == maxdist)
    
    init_local_ctrl_pts_pos = np.array([indx[0][0],indx[1][0]]) # take only the first element (two pts)
    init_local_ctrl_pts = vecmin.vecmin(np.reshape(init_local_ctrl_pts_pos.T,[1,-1])[0])
    init_local_ctrl_pts_dist_table = ctrl_pts2dist_table(init_local_ctrl_pts,local_dist_table)
    
    ## Search & Check routine
    local_ctrl_pts_num = 2
    local_ctrl_pts = init_local_ctrl_pts.copy()
    local_ctrl_pts_dist_table = init_local_ctrl_pts_dist_table.copy()
    global_dist_table_half = dist_table_shriek(global_dist_table)
    
    hitflag = 0 # Unique candidate found
    failflag = 0 # Search/check failed
    
    while hitflag == 0 and failflag == 0:
        
        search_candidates,search_candidates_prev = dist_table_check(local_ctrl_pts_dist_table,global_dist_table,T)
        
        hit = searchcand()
        candnum = len(search_candidates)
        if candnum == 0: # No candidate found
            failflag = 1 # Search failed
            break
        elif candnum == 1:
            hitflag = 1 # Search found
            hit = search_candidates[0]
            break
        # elif local_ctrl_pts_num == local_CtrlPtNum: # All local CtrlPts known has been used
        #     hitflag = 0
        #     
        #     semihitflag = 1
        #     # eliminate the identical candidates case
        #     for ii in range(candnum):
        #         cand_ctrl_pts_vecmin_a = search_candidates[ii].ctrl_pts_minvec
        #         for jj in range(ii+1,candnum):
        #             cand_ctrl_pts_vecmin_b = search_candidates[jj].ctrl_pts_minvec
        #             if np.nansum(np.abs(cand_ctrl_pts_vecmin_a-cand_ctrl_pts_vecmin_b)) == 0: # Same points
        #                 continue
        #             else:
        #                 semihitflag = 0 # Search failed
        #                 break
        #         if semihitflag = 0:
        #             break 
        #     # Further expansion of this section might lead to better handling of special cases (small detected CtrlPts number)
                    
        # Add a new control point from local
        if hitflag == 0 or failflag == 0: # Search to be continued
            local_ctrl_pts_num += 1
            local_ctrl_pts = add_local_ctrl_points(local_ctrl_pts,local_CtrlPtNum)
            local_ctrl_pts_dist_table = ctrl_pts2dist_table(local_ctrl_pts,local_dist_table)
        
    return hit,local_ctrl_pts,search_candidates,search_candidates_prev
    
def mutate_array(a,ii,jj,axis=0):
    if axis == 0: # col
        tmp = a[:,ii].copy()
        a[:,ii] = a[:,jj]
        a[:,jj] = tmp.copy()
    elif axis == 1: # line
        tmp = a[ii,:].copy()
        a[ii,:] = a[jj,:]
        a[jj,:] = tmp.copy()
    
    return a

def swap_vec(v,ii,jj):
    v_swap = v.copy()
    tmp = v_swap[ii].copy()
    v_swap[ii] = v_swap[jj]
    v_swap[jj] = tmp.copy()
    
    return v_swap
    
def candidate_center(search_candidate,CtrlPt_pos):
    
    ctrl_pts = search_candidate.ctrl_pts
    x = CtrlPt_pos[0,ctrl_pts]
    y = CtrlPt_pos[1,ctrl_pts]
    
    candcenter = np.zeros(2)
    candcenter[0] = np.mean(x)
    candcenter[1] = np.mean(y)
    
    return candcenter

def prev_pos_check(search_candidates,CtrlPt_pos,prev_pos,prev_pos_TH):
    
    search_candidates_filtered = []
    for ii in range(len(search_candidates)):
        candcenter = candidate_center(search_candidates[ii],CtrlPt_pos)
        if np.sum(np.abs(candcenter-prev_pos))<prev_pos_TH:
            search_candidates_filtered.append(search_candidates[ii])
    
    return search_candidates_filtered

## Control Points Template Match with outliner handling
def control_points_template_match_outH(local_dist_table,global_dist_table,minpts=4,maxpts=7,T=0.2):
# The algorithm assumes that at least minpts real CtrlPts are detected
# The algorithm assumes that if trustpts CtrlPts are detected, they are all real
# Given that no trustpts CtrlPts pass the check, we further assume that trustpts -= 1 CtrlPts are trustworthy
# In the case of multiple candidates passing search and check (different combinations), return the first one found

    # randomly select all to minpts for match
    local_CtrlPt_num = len(local_dist_table[0,:])
    trustpts = np.min([local_CtrlPt_num,maxpts]) # We trust trustpts = maxpts points to be robust enough
    
    ctrl_pts = []
    hit = []
    search_candidates = []
    
    while trustpts >= minpts:
        
        all_combinations = comb_index(local_CtrlPt_num,trustpts) # Create all mutation
        combination_outliner_possibility = np.zeros(len(all_combinations))
        
        hitflag = 0 # Hit, unique candidate
        multihitflag = 0 # multiple candidates detected
        
        while hitflag == 0 and multihitflag == 0:
            if np.nanmin(combination_outliner_possibility) == np.inf:
                break # All combinations checked
                
            indx = np.nanargmin(combination_outliner_possibility) # Find the combination that is most unlikely to contain outliners
            combination_sel = all_combinations[indx,:].copy() # get the combination
            
            # Get the local_dist_table_sel from local_dist_table
            combination_sel_dist_table = ctrl_pts2dist_table(combination_sel,local_dist_table)
            
            # Rearrange the combination: max distance first search
            maxdist = np.nanmax(combination_sel_dist_table)
            combination_sel_dist_table_half = dist_table_shriek(combination_sel_dist_table)
            indx = np.where(combination_sel_dist_table_half == maxdist)
            init_ptA = indx[0][0]
            init_ptB = indx[1][0] # Get the points indx
            
            local_CtrlPt_sel = swap_vec(combination_sel,0,init_ptA)
            local_CtrlPt_sel = swap_vec(local_CtrlPt_sel,1,init_ptB)
            local_dist_table_sel = ctrl_pts2dist_table(local_CtrlPt_sel,local_dist_table)
            
            # Find/check if the candidates could be found in the global map
            search_candidates,search_candidates_prev = dist_table_check(local_dist_table_sel,global_dist_table,T)
            
            candnum = len(search_candidates)
            if candnum == 1:
                search_candidates[0].set_local_pts(local_CtrlPt_sel)
                hit = search_candidates[0]
                hitflag = 1
                break
            elif candnum>1: # More than one candidates found
                for ii in range(candnum):
                    search_candidates[ii].set_local_pts(local_CtrlPt_sel)
                multihitflag = 1 # Multi-hit occurs.Lower level search no longer necessary, as the algorithm has produced multi-ple equally likely candidates
                break
        
            # Update the combination_outliner_possibility
            for ii in range(trustpts):
                outliner_possible = combination_sel[ii]
                for jj in range(trustpts):
                    combination_dim_vec = all_combinations[:,jj]
                    indx = np.where(combination_dim_vec == outliner_possible)
                    combination_outliner_possibility[indx] = combination_outliner_possibility[indx] + 1
                    
            # Find the perfect match, assign np.inf to its outliner possibility
            indx = np.where(np.nansum(np.abs(all_combinations - combination_sel),axis=1) == 0)
            combination_outliner_possibility[indx] = np.inf
            
        if hitflag == 0 and multihitflag == 0: # No hit or multi-hit
            trustpts = trustpts - 1 # reduce to the next level
        else:
            break
    
    if hitflag == 1 and trustpts == maxpts: # In the case of hit with maxpts, expand
        
        local_CtrlPt_sel_l = local_CtrlPt_sel.tolist()
        
        # Expand local_CtrlPt_sel to the highest possible number
        # Assume that given the high number of inliners, any outliner would be eliminated and any inliner checked in would be fine.
        for ii in range(local_CtrlPt_num):
            indx = np.where(local_CtrlPt_sel-ii==0)
            if not indx: # Not empty
                local_CtrlPt_sel_l.append(ii)
            else:
                continue
            
            local_CtrlPt_sel_trial = np.array(local_CtrlPt_sel_l)
            local_dist_table_sel = ctrl_pts2dist_table(local_CtrlPt_sel_trial,local_dist_table)
            search_candidates,search_candidates_prev = dist_table_check(local_dist_table_sel,global_dist_table,T)
            
            candnum = len(search_candidates)
            if candnum == 1:
                hit = search_candidates[0]
                local_CtrlPt_sel = local_CtrlPt_sel_trial.copy() # Update local CtrlPt selection
            
    return hit,local_CtrlPt_sel,search_candidates,search_candidates_prev

## With angle check and outliner handling
def control_points_template_match(local_dist_table,local_ang_table,CtrlPt_pos,global_dist_table,global_ang_table,minpts=4,maxpts=7,disttol=0.2,angtol=5/180*pi,prev_pos=[],prev_pos_TH= 1):
# Angle check version
# The algorithm assumes that at least minpts real CtrlPts are detected
# The algorithm assumes that if trustpts CtrlPts are detected, they are all real
# Given that no trustpts CtrlPts pass the check, we further assume that trustpts -= 1 CtrlPts are trustworthy
# In the case of multiple candidates passing search and check (different combinations), return the first one found

    # randomly select all to minpts for match
    local_CtrlPt_num = len(local_dist_table[0,:])
    trustpts = np.min([local_CtrlPt_num,maxpts]) # We trust trustpts = maxpts points to be robust enough
    
    ctrl_pts = []
    hit = []
    search_candidates = []
    
    while trustpts >= minpts:
        
        all_combinations = comb_index(local_CtrlPt_num,trustpts) # Create all mutation
        combination_outliner_possibility = np.zeros(len(all_combinations))
        
        hitflag = 0 # Hit, unique candidate
        multihitflag = 0 # multiple candidates detected
        
        while hitflag == 0 and multihitflag == 0:
            if np.nanmin(combination_outliner_possibility) == np.inf:
                break # All combinations checked
                
            indx = np.nanargmin(combination_outliner_possibility) # Find the combination that is most unlikely to contain outliners
            combination_sel = all_combinations[indx,:].copy() # get the combination
            
            # Get the local_dist_table_sel from local_dist_table
            combination_sel_dist_table = ctrl_pts2dist_table(combination_sel,local_dist_table)
            
            # Rearrange the combination: max distance first search
            maxdist = np.nanmax(combination_sel_dist_table)
            combination_sel_dist_table_half = dist_table_shriek(combination_sel_dist_table)
            indx = np.where(combination_sel_dist_table_half == maxdist)
            init_ptA = indx[0][0]
            init_ptB = indx[1][0] # Get the points indx
            
            local_CtrlPt_sel = swap_vec(combination_sel,0,init_ptA)
            local_CtrlPt_sel = swap_vec(local_CtrlPt_sel,1,init_ptB)
            local_dist_table_sel = ctrl_pts2dist_table(local_CtrlPt_sel,local_dist_table)
            
            # Find/check if the candidates could be found in the global map
            search_candidates,search_candidates_prev = dist_table_check(local_dist_table_sel,global_dist_table,disttol)
            
            ## Angle check
            candnum = len(search_candidates)
            if candnum >= 1:
                for ii in range(candnum):
                    search_candidates[ii].set_local_pts(local_CtrlPt_sel)
                search_candidates = ang_table_check(search_candidates,global_ang_table,local_ang_table,angtol)
            
            ## pos_previous filtering
            candnum = len(search_candidates)
            if candnum >= 1:
                if not prev_pos:
                    pass
                else:
                    search_candidates = prev_pos_check(search_candidates,CtrlPt_pos,prev_pos,prev_pos_TH)
                
            candnum = len(search_candidates)
            if candnum == 1:
                hit = search_candidates[0]
                hitflag = 1
                break
            elif candnum>1: # More than one candidates found
                multihitflag = 1 # Multi-hit occurs.Lower level search no longer necessary, as the algorithm has produced multi-ple equally likely candidates
                break
        
            # Update the combination_outliner_possibility
            for ii in range(trustpts):
                outliner_possible = combination_sel[ii]
                for jj in range(trustpts):
                    combination_dim_vec = all_combinations[:,jj]
                    indx = np.where(combination_dim_vec == outliner_possible)
                    combination_outliner_possibility[indx] = combination_outliner_possibility[indx] + 1
                    
            # Find the perfect match, assign np.inf to its outliner possibility
            indx = np.where(np.nansum(np.abs(all_combinations - combination_sel),axis=1) == 0)
            combination_outliner_possibility[indx] = np.inf
            
        if hitflag == 0 and multihitflag == 0: # No hit or multi-hit
            trustpts = trustpts - 1 # reduce to the next level
        else:
            break
    
    if hitflag == 1 and trustpts == maxpts: # In the case of hit with maxpts, expand
        
        local_CtrlPt_sel_l = local_CtrlPt_sel.tolist()
        
        # Expand local_CtrlPt_sel to the highest possible number
        # Assume that given the high number of inliners, any outliner would be eliminated and any inliner checked in would be fine.
        for ii in range(local_CtrlPt_num):
            indx = np.where(local_CtrlPt_sel-ii==0)
            if not indx: # Not empty
                local_CtrlPt_sel_l.append(ii)
            else:
                continue
            
            local_CtrlPt_sel_trial = np.array(local_CtrlPt_sel_l)
            local_dist_table_sel = ctrl_pts2dist_table(local_CtrlPt_sel_trial,local_dist_table)
            search_candidates,search_candidates_prev = dist_table_check(local_dist_table_sel,global_dist_table,T)
            
            candnum = len(search_candidates)
            if candnum == 1:
                hit = search_candidates[0]
                local_CtrlPt_sel = local_CtrlPt_sel_trial.copy() # Update local CtrlPt selection

    return hit,local_CtrlPt_sel,search_candidates

