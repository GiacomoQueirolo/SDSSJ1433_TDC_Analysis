#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pickle as pkl
import numpy as np
from pycs3.tdcomb.comb import Group
import numpy as np
import matplotlib.pyplot as plt
from corner import quantile
import copy
import warnings


# In[3]:


class Error():
    
    def __init__(self,error_path=None):
        self.error_path = error_path
        if os.path.exists(error_path):
            self.create_error()
        
    def get_dt(self,meas_or_sim):
        if self.error_path is None:
            raise RuntimeError("Give me a path to the simulations results")
        list_file = os.listdir(path=self.error_path)
        list_file = [list_file[i] for i in range(len(list_file)) if "runresults.pkl" in list_file[i]]
        if len(list_file)==0:
            raise RuntimeError("No resulting files in this directory: "+str(self.error_path))
        for i,file in enumerate(list_file):
            rr = pkl.load(open(self.error_path+"/"+file,"rb"))
            if meas_or_sim=="meas":
                tsarr = rr.tsarray
            elif meas_or_sim=="sim":
                tsarr = rr.truetsarray
            # MOD_DELTA
            # do not consider tA,tB and tC, but DtAB,DtAC AND DtBC 
            # tsarr.shape = 20,3
            tsarrT = tsarr.T # shape 3,20
            new_tsarrT = [tsarrT[1:]-tsarrT[0]] # AB,AC: shape 2,20
            new_tsarrT.append(tsarr[2]-tsarr[1]) # BC: shape 3,20
            new_tsarr = np.transpose(new_tsarr) # shape 20,3 
            if i==0:
                dts = new_tsarr
            else:
                dts = np.vstack((dts,new_tsarr))
        return np.transpose(dts) #shape: dim, n*steps

    def get_meas_dt(self):
        if hasattr(self,"measdt") is False:
            self.measdt = self.get_dt("meas")
        self.labels_measdt = ["AB","AC","BC"]
        return self.measdt
    
    def get_sim_dt(self):
        if hasattr(self,"simdt") is False:
            self.simdt = self.get_dt("sim")
        self.labels_simdt = ["AB","AC","BC"]
        return self.simdt 

    def get_distr(self):
        if hasattr(self,"distr"):
            return self.distr
        elif not hasattr(self,"measdt") and not hasattr(self,"simdt"):
            self.meas = self.get_meas_dt()
            self.sim  = self.get_sim_dt()
        #print("Note:Error is defined as measured dt - simulated dt")
        distr = self.measdt - self.simdt
        self.distr = distr
        if self.labels_measdt==self.labels_simdt:
            self.labels = self.labels_simdt
        else:
            raise RuntimeError("Check the labels/dimensions")
        return self.distr
            
    def create_error(self):    
        err_all = self.get_distr()

        errors_up,errors_down,rnd,sys =[],[],[],[]
        tot=[]
        for i in range(len(err_all)):
            min_e,med_e,max_e = quantile(err_all[i],q=[0.16,0.5,0.84])
            err_up,err_down = max_e-med_e,med_e-min_e
            errors_up.append(err_up)
            errors_down.append(err_down)
            rnd.append((err_up+err_down)/2.)
            sys.append(med_e)
        #tot =  sqrt(rnd² +sys²)
        ############################################################
        #tot = [np.sqrt(rnd[i]**2 + sys[i]**2) for i in range(len(err_all))]
        tot = sqrt_sum(rnd,sys)
        ############################################################
        self.err_up   = errors_up
        self.err_down = errors_down
        self.rnd = rnd
        self.sys = sys
        self.tot = tot
        self.accuracy = np.mean(tot) #prob not needed anymore


# In[2]:


class Group():
    def __init__(self,error=None,data_distr=None,labels=["AB","AC","BC"],name="", color="royalblue"):
        
        self.data   = data_distr
        self.labels = labels
        self.name   = name
        self.color  = color
        
        if self.data is not None:
            self.results = np.mean(data_distr,axis=1)  # as usual median or mean? 
        else:
            self.results = None
        # All the error description are inherited by the Errorbar class
        if error is not None:
            if isinstance(error,Error):
                self.error = error
            elif type(error)==str : 
                error       = Error(error)
                self.error  = error
            else:
                error_tmp       = Error("None")
                error_tmp.distr = error 
                self.error      = error_tmp
            self.error.create_error()
            # check that we have the right dataset:
            if labels != self.error.labels:
                raise RuntimeError("Check the labels/dimensions!")
            self.err_distr  = self.error.get_distr()
            self.err_up     = self.error.err_up
            self.err_down   = self.error.err_down
            self.rnd_error  = self.error.rnd
            self.sys_error  = self.error.sys
            self.tot_error  = self.error.tot
            self.accuracy   = self.error.accuracy


# In[ ]:


## Useful functions 
        
def getresults(data_path,name="",error=None,color="royalblue",labels=["AB","AC","BC"]):
    """
    ;My version of pycs3.tdcomb.comb.getresults; no need for CScontainer
    """

    with open(str(data_path)+"/td.data","rb") as f:
        timedelays = pkl.load(f)
    ###############################
    #To better correct one day
    if timedelays[-1]==0:
        timedelays=timedelays[:-1]
    ###############################
    timedelays = np.transpose(timedelays) #shape now: lcs, mci
    # create a Group out of them
    group =  Group(data_distr=timedelays,error = error,name=name,labels=labels, color=color)
    return group

def sqrt_sum(a,b):
    if np.shape(a)==():
        return np.sqrt(a**2 + b**2)
    else:
        if len(a)!=len(b):
            raise ValueError("Len of a and b must be equal, not: "+str(len(a))+" and "+str(len(b)))
        return [np.sqrt(a[i]**2 + b[i]**2) for i in range(len(a))]
    
def combine_groups(G1,G2):
    if G1.labels!=G2.labels: #pragma: no cover
        raise RuntimeError("Give me two groups of the same dataset!")
    
    #Following method A in Notes - 27th October 
    #res_comb = np.mean(np.vstack([G1.data,G2.data]),axis=0)                                      
    res_comb = np.mean([G1.results,G2.results],axis=0)                                          
    err_comb = sqrt_sum(G1.tot_error,G2.tot_error)
    
    #comb_error       = Error("none")
    #comb_error.distr = None
    #comb_error.tot   = err_comb
    
    G_comb           = Group(error=None,labels=G1.labels)
    G_comb.results   = res_comb
    G_comb.tot_error = err_comb
    
    G_comb.color = "black"
    G_comb.name  = "Combined "+G1.name.replace("Combined ","")+" and "+G2.name.replace("Combined ","")
    return G_comb

##########################
def get_ref_index(series):
    err = [G.tot_error for G in series]
    ref_index = int(np.where(err==np.min(err))[0])
    return ref_index

def tau(dt_i,dt_j,sig_i,sig_j): 
    #simplified version bc considering method A
    tau_val = abs(dt_i-dt_j)/sqrt_sum(sig_i,sig_j) #~ Z val
    return tau_val

def tau_G(G1,G2):
    if G1.labels!=G2.labels: #pragma: no cover
        raise RuntimeError("Give me two groups of the same dataset!")
    tau_G = []
    for i in range(len(G1.labels)):
        dt1  = G1.results[i]
        dt2  = G2.results[i]
        sig1 = G1.tot_error[i]
        sig2 = G2.tot_error[i]
        tau_G.append(tau(dt1,dt2,sig1,sig2))
    return max(tau_G)

def combine_series(series,sigmathresh=0.5,return_combined_list=False):
    series = copy.copy(series)
    # find the reference group (most accurate overall)
    ref_index = get_ref_index(series)
    ref_G     = copy.copy(series[ref_index])
    ignore_G  = [ref_G.name]
    print("Initial best result: ",ref_G.name)
    # check tension
    tension_series =[] 
    for G in series:
        if G.name not in ignore_G and tau_G(G,ref_G)>=sigmathresh:
            tension_series.append(G)
    
    print("Combining results series...")
    while len(tension_series)>0:
        tns_index = get_ref_index(tension_series)
        tns_ref   = tension_series[tns_index]
        ignore_G.append(tns_ref.name)
        ref_G = combine_groups(ref_G,tns_ref)
        tension_series =[]
        for G in series:
            if G.name not in ignore_G and tau_G(G,ref_G)>=sigmathresh:
                tension_series.append(G)
    print("Combined:",ignore_G)
    ref_G.name=r"Combined result with $\tau_{thresh}=$"+str(sigmathresh)
    ref_G.combined_names = ignore_G
    if return_combined_list:
        return ref_G,ignore_G
    else:
        return ref_G


# In[ ]:



### Try: method B (see Nov 2nd-3rd)
def corr_distr(err_distr,dt_meas):
    if len(np.transpose(err_distr))==len(dt_meas):
        err_distr=np.transpose(err_distr)
    elif not len(err_distr)==len(dt_meas):
        raise RuntimeError("err_distr and dt_meas must have at least 1 dim in common,                            instead ",np.shape(err_distr)," and ",np.shape(dt_meas))
    corr_dist,sys = [],[]
    for i in range(len(dt_meas)):
        sys.append(np.mean(err_distr[i]))
        corr_dist.append(err_distr[i]-sys[-1]+dt_meas[i])
    return np.array(corr_dist),np.array(np.abs(sys))

    
def combine_groups_methodB(G1,G2):
    if G1.labels!=G2.labels: #pragma: no cover
        raise RuntimeError("Give me two groups of the same dataset!")
    
    #Following method B in Notes - 27th October 
    
    #first we "correct" the error distr
    G1.corr_distr,G1.sys = corr_distr(G1.err_distr,G1.results)
    G2.corr_distr,G2.sys = corr_distr(G2.err_distr,G2.results)
    
    #we now can combine them and obtain the dt resulting from the combine distr
    #and compute the combined sys
    comb_distr = np.hstack([G1.corr_distr,G2.corr_distr])
    comb_res   = np.mean(comb_distr,axis=1) # always to consider if median or mean
    comb_sys   = sqrt_sum(G1.sys,G2.sys)
    
    #now we restore the error distr. as centered on the sys,
    #while remembering the time result
    comb_error       = Error()
    comb_error.distr = [comb_distr[k] - comb_res[k] + comb_sys[k] for k in range(len(comb_sys))]
    G_comb           = Group(data_distr = comb_distr,                             error  = comb_error,                             labels = G1.labels)
    
    G_comb.color = "black"
    G_comb.name  = r"Combined "+G1.name.replace("Combined ","")+" and "+G2.name.replace("Combined ","")
    return G_comb

def get_sig_prime(group,i,prime=""):
    # see Notes 3rd Nov. '21
    sig_up  = group.err_up[i]
    sig_low = group.err_down[i]
    sig_sys = group.sys_error[i]
    if prime=="up":
        return np.sqrt(sig_up**2 + sig_low*sig_up + 2*sig_sys**2)
    elif prime=="low":
        return np.sqrt(sig_low**2 + sig_low*sig_up + 2*sig_sys**2)
    else:
        raise ValueError("Either upper or lower error")
                
def tau_G_methodB(G1,G2):
    if G1.labels!=G2.labels: #pragma: no cover
        raise RuntimeError("Give me two groups of the same dataset!")
    tau_G = []
    for i in range(len(G1.labels)):
        dt1  = G1.results[i]
        dt2  = G2.results[i]
        if dt1>dt2:
            sig1 = get_sig_prime(G1,i,prime="low")
            sig2 = get_sig_prime(G2,i,prime="up")
        else:
            sig1 = get_sig_prime(G1,i,prime="up")
            sig2 = get_sig_prime(G2,i,prime="low")
        tau = abs(dt1-dt2)/sqrt_sum(sig1,sig2)
        tau_G.append(tau)
    return max(tau_G)

def combine_series_methodB(series,sigmathresh=0.5,return_combined_list=False):
    warnings.warn("\nWARNING: we are here using the method B to combine groups.     See notes of 27thOct.\n",)
    series = copy.copy(series)
    #Set all the group in the series to have the same bin
    #comb_bins = get_bins(series,n_bins=150,return_bins=True) 
    #set all the series' group such that they are shifted and have their (shifted) density distribution
    # find the reference group (most accurate overall)
    ref_index = get_ref_index(series) # this, provided that each Group has its own tot_error,doens't change
    ref_G     = copy.copy(series[ref_index])
    ignore_G  = [ref_G.name]
    print("Initial best result: ",ref_G.name)
    # check tension
    tension_series =[] 
    for G in series:
        if G.name not in ignore_G and tau_G_methodB(G,ref_G)>=sigmathresh:
            tension_series.append(G)
    
    print("Combining results series...")
    while len(tension_series)>0:
        tns_index = get_ref_index(tension_series)
        tns_ref   = tension_series[tns_index]
        ignore_G.append(tns_ref.name)
        ref_G = combine_groups_methodB(ref_G,tns_ref)
        tension_series =[]
        for G in series:
            if G.name not in ignore_G and tau_G_methodB(G,ref_G)>=sigmathresh:
                tension_series.append(G)
    print("Combined:",ignore_G)
    ref_G.name+=r"\nCombined result with $\tau_{thresh}=$"+str(sigmathresh)
    ref_G.combined_names = ignore_G
    if return_combined_list:
        return ref_G,ignore_G
    else:
        return ref_G

