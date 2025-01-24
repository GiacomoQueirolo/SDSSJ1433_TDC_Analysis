#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import copy
import warnings
import numpy as np
import pickle as pkl
from corner import quantile
from pycs3.tdcomb.comb import Group

from Utils.tools import get_analdir
from Utils.math_tools import sqrt_sum, sqrt_sum_list

class Error_main():
    
    def __init__(self,error_path=None):
        self.error_path = error_path
        if error_path is not None:
            if os.path.exists(error_path):
                self.create_error()   
    def get_distr(self):
        if hasattr(self,"distr"):
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
        
class Error():
    def __init__(self,error_path=None,wD=False):
        self.error_path = error_path
        self.wD         = wD
        if error_path is not None:
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
            if not self.wD:
                new_tsarrT = np.array([tsarrT[1]-tsarrT[0],tsarrT[2]-tsarrT[0],tsarrT[2]-tsarrT[1]]) # AB,AC,BC: shape 3,20
            else:
                new_tsarrT = np.array([tsarrT[1]-tsarrT[0],tsarrT[2]-tsarrT[0],tsarrT[3]-tsarrT[0]]) # AB,AC,AD: shape 3,20
            new_tsarr = np.transpose(new_tsarrT) # shape 20,3 
            if i==0:
                dts = new_tsarr
            else:
                dts = np.vstack((dts,new_tsarr))
        return np.transpose(dts) #shape: dim, n*steps

    def get_meas_dt(self):
        if hasattr(self,"measdt") is False:
            self.measdt = self.get_dt("meas")
        if self.wD:
            self.labels_measdt = ["AB","AC","AD"]
        else:
            self.labels_measdt = ["AB","AC","BC"]
        return self.measdt
    
    def get_sim_dt(self):
        if hasattr(self,"simdt") is False:
            self.simdt = self.get_dt("sim")
        if self.wD:
            self.labels_simdt = ["AB","AC","AD"]
        else:
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
            print(self.labels_measdt,self.labels_simdt)
            raise RuntimeError("Check the labels/dimensions")
        return self.distr
            
    def create_error(self):    
        err_all = self.get_distr()

        errors_up,errors_down,rnd,sys,tot =[],[],[],[],[]
        # TEST : MOD_INTRSCATTER 
        data_path = get_analdir(self.error_path)
        if "/"==data_path[-1]:
            data_path=data_path[:-1]
        data_path = "/".join(data_path.split("/")[:-1])
        if str(data_path)=="":
            intrsct = [0,0,0]
        else:
        
            try:
                Gr = getresults(data_path=data_path)
            except:
                print(self.error_path)
                exit()
            intrsct = np.std(Gr.data,axis=1)
        print("intrinsic error",intrsct)        
         
        for i in range(len(err_all)):
            min_e,med_e,max_e = quantile(err_all[i],q=[0.16,0.5,0.84])
            err_up,err_down = max_e-med_e,med_e-min_e
            errors_up.append(err_up)
            errors_down.append(err_down)
            rnd.append((err_up+err_down)/2.)
            sys.append(np.sign(med_e)*sqrt_sum(med_e,intrsct[i]))
            #print("added to sys:",med_e,"->",sys[i])
            #sys.append(med_e)
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


class Error_mag():
    
    def __init__(self,error_path=None,wD=False):
        self.error_path = error_path
        self.wD         = wD
        if error_path is not None:
            if os.path.exists(error_path):
                self.create_error()
        
    def get_dmag(self,meas_or_sim):
        if self.error_path is None:
            raise RuntimeError("Give me a path to the simulations results")
        list_file = os.listdir(path=self.error_path)
        list_file = [list_file[i] for i in range(len(list_file)) if "runresults.pkl" in list_file[i]]
        if len(list_file)==0:
            raise RuntimeError("No resulting files in this directory: "+str(self.error_path))
        for i,file in enumerate(list_file):
            rr = pkl.load(open(self.error_path+"/"+file,"rb"))
            if meas_or_sim=="meas":
                magsarr = rr.magsarray
                if self.wD:
                    self.labels_measdmag  = ["AB","AC","AD"]
                else:
                    self.labels_measdmag = ["AB","AC","BC"]
            elif meas_or_sim=="sim":
                magsarr = rr.truemagsarray
                if self.wD:
                    self.labels_simdmag  = ["AB","AC","AD"]
                else:
                    self.labels_simdmag = ["AB","AC","BC"]
            elif meas_or_sim=="error":
                # error defined as dmag- simulated dmag
                magsarr = rr.magsarray - rr.truemagsarray
                if self.wD :
                    self.labels          = ["AB","AC","AD"]
                    self.labels_simdmag  = ["AB","AC","AD"]
                    self.labels_measdmag = ["AB","AC","AD"]
                else:
                    self.labels          = ["AB","AC","BC"]
                    self.labels_simdmag  = ["AB","AC","BC"]
                    self.labels_measdmag = ["AB","AC","BC"]
            # MOD_DELTA
            # do not consider mA,mB and mC, but DmAB,DmAC AND DmBC 
            # magsarr.shape = 20,3
            magsarrT = magsarr.T # shape 3,20
            if not self.wD :
                Dmagsarr = np.array([magsarrT[1]-magsarrT[0],magsarrT[2]-magsarrT[0],magsarrT[2]-magsarrT[1]]) # AB,AC,BC: shape 3,20
            else:
                Dmagsarr = np.array([magsarrT[1]-magsarrT[0],magsarrT[2]-magsarrT[0],magsarrT[3]-magsarrT[0]]) # AB,AC,AD
            #new_magsarr = np.transpose(new_magsarrT) # shape 20,3 
            if i==0:
                dms = Dmagsarr
            else:
                dms = np.hstack((dms,Dmagsarr))
        return dms #shape: dim, n*steps

    def get_meas_dmag(self):
        if hasattr(self,"measdmag") is False:
            self.measdmag = self.get_dmag("meas")
        return self.measdmag
    
    def get_sim_dmag(self):
        if hasattr(self,"simdmag") is False:
            self.simdmag = self.get_dmag("sim")
        return self.simdmag 

    def get_distr(self):
        if hasattr(self,"distr"):
            return self.distr
        #print("Note:Error is defined as measured dmag - simulated dmag")
        distr = self.get_dmag("error")#self.get_meas_dmag() - self.get_sim_dmag()
        self.distr = distr
        if self.labels_measdmag==self.labels_simdmag:
            self.labels = self.labels_simdmag
        else:
            print("mag",self.labels_measdmag,self.labels_simdmag)
            raise RuntimeError("Check the labels/dimensions")
        return self.distr
            
    def create_error(self):    
        err_all = self.get_distr()

        errors_up,errors_down,rnd,sys =[],[],[],[]
        # TEST : MOD_INTRSCATTER 
        data_path = get_analdir(self.error_path)
        if "/"==data_path[-1]:
            data_path=data_path[:-1]
        data_path = "/".join(data_path.split("/")[:-1])
        if str(data_path)=="":
            # measn that we are doing combined error,
            # intrinsic scatter already considered
            intrsct = [0,0,0]
        else:
            Gr = getresults_mag(data_path=data_path)
            intrsct = np.std(Gr.data,axis=1)
        print("intrinsic error mag",intrsct)        
    
        tot=[]
        for i in range(len(err_all)):
            min_e,med_e,max_e = quantile(err_all[i],q=[0.16,0.5,0.84])
            err_up,err_down = max_e-med_e,med_e-min_e
            errors_up.append(err_up)
            errors_down.append(err_down)
            rnd.append((err_up+err_down)/2.)
            sys.append(np.sign(med_e)*sqrt_sum(med_e,intrsct[i]))
            #sys.append(med_e)
            #print("added to sys:",med_e,"->",sys[i])
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



class Group():
    def __init__(self,error=None,data_distr=None,labels=["AB","AC","BC"],name="", color="royalblue",\
                 config=None,knt=None,mltype=None,ml_config=None):
        
        self.data      = data_distr
        self.labels    = labels
        self.name      = name
        self.color     = color
        self.wD        = False if "BC" in self.labels else True
        self.config    = config
        self.knt       = knt
        self.mltype    = mltype
        self.ml_config = ml_config

        if self.data is not None:
            self.results = np.mean(data_distr,axis=1)  # as usual median or mean? 
        else:
            self.results = None
        # All the error description are inherited by the Errorbar class
        if error is not None:
            if isinstance(error,(Error,Error_mag)):
                self.error = error
            elif type(error)==str : 
                error       = Error(error,wD=self.wD)
                self.error  = error
            else:
                error_tmp       = Error("None",wD=self.wD)
                error_tmp.distr = error 
                self.error      = error_tmp
            self.error.create_error()
            # check that we have the right dataset:
            if labels != getattr(self.error,"labels",labels):
                raise RuntimeError("Check the labels/dimensions!")
            self.err_distr  = self.error.get_distr()
            self.err_up     = self.error.err_up
            self.err_down   = self.error.err_down
            self.rnd_error  = self.error.rnd
            self.sys_error  = self.error.sys
            self.tot_error  = self.error.tot
            self.accuracy   = self.error.accuracy





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

     
def getresults_mag(data_path,name="",error=None,color="royalblue",labels=["AB","AC","BC"]):
    """
    ;same as getresults, but for mag
    """
    if "splml" in str(data_path):
        print("Not implemented for spline ML")
        return None
    with open(str(data_path)+"/dmags.data","rb") as f:
        dmags = pkl.load(f)
    if labels[-1]=="BC":
        # to implement better:
        dmags = np.transpose(dmags) #shape now: lcs, mci
        # create a Group out of them
        DmagBC = dmags[2] - dmags[1] # BC
        Dmags  = dmags[1:]- dmags[0] #rel to A : AB,AC
        Dmags  = np.append(Dmags,[DmagBC],0).tolist()
        # now AB,AC,BC -> mostly to be consistent w. dt
    else:
        Dmags = np.transpose(dmags).tolist()
    group =  Group(data_distr=Dmags,error = error,name=name,labels=labels, color=color)
    return group


    
###############################
# archival function, don't use
def combine_group_list(list_G):
    raise RuntimeWarning("Don't use this one, go for the methodB")
    for Gi in list_G[1:]:
        if Gi.labels!=list_G[0].labels: #pragma: no cover
            raise RuntimeError("Give me two groups of the same dataset!")
    #Following method A in Notes - 27th October 
    res_comb = np.mean([Gi.results for Gi in list_G],axis=0)                                          
    err_comb = sqrt_sum_list([Gi.tot_error for Gi in list_G])/len(list_G)
        
    G_comb           = Group(error=None,labels=G1.labels)
    G_comb.results   = res_comb
    G_comb.tot_error = err_comb
    
    G_comb.color = "black"
    G_comb.name  = "Combined "+str([Gi.name.replace("Combined ","")+" and " for Gi in list_G])
    return G_comb
###############################

def combine_group_list_methodB(list_G):
    if len(list_G)==1:
        print("Only one element in the given list to combine")
        return list_G[0]
    for Gi in list_G[1:]:
        if Gi.labels!=list_G[0].labels: #pragma: no cover
            raise RuntimeError("Give me two groups of the same dataset!")
        G1=list_G[0]
    #Following method B in Notes - 27th October 
    
    #first we "correct" the error distr
    for Gi in list_G:
        Gi.corr_distr,Gi.sys = corr_distr(Gi.err_distr,Gi.results)
    comb_distr = np.hstack([Gi.corr_distr for Gi in list_G])
    comb_res   = np.mean(comb_distr,axis=1) # always to consider if median or mean
    comb_sys   = sqrt_sum_list([Gi.sys for Gi in list_G])
    
    #now we restore the error distr. as centered on the sys,
    #while remembering the time result
    comb_error       = Error()
    comb_error.distr = [comb_distr[k] - comb_res[k] + comb_sys[k] for k in range(len(comb_sys))]
    G_comb           = Group(data_distr = comb_distr, error = comb_error, labels = G1.labels)

    G_comb.color = "black"
    G_comb.name  = "Combined "+str([Gi.name.replace("Combined ","")+" and " for Gi in list_G])
    return G_comb


def get_ref_index(series):
    #err = [G.tot_error for G in series]
    # 16th dec 22 : this way this is selecting only the result with the single lower error
    # while we might want the one with - in average - the lower error
    # shouldn't change too much
    err       = [np.mean(G.tot_error) for G in series]
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

#######################################################################
# archival function. don't use
def combine_series(series,sigmathresh=0.5,return_combined_list=False):
    raise RuntimeWarning("Archival function, don't use. Use methodB")
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
#######################################################################


### Try: method B (see Nov 2nd-3rd)
def corr_distr(err_distr,dt_meas):
    if len(np.transpose(err_distr))==len(dt_meas):
        err_distr=np.transpose(err_distr)
    elif not len(err_distr)==len(dt_meas):
        raise RuntimeError("err_distr and dt_meas must have at least 1 dim in common, instead ",np.shape(err_distr)," and ",np.shape(dt_meas))
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
    G_comb           = Group(data_distr = comb_distr, error = comb_error, labels = G1.labels)
    
    G_comb.color = "black"
    G_comb.name  = r"Combined "+G1.name.replace("Combined ","")+" and "+G2.name.replace("Combined ","")
    return G_comb

def get_sig_prime(group,i,prime=""):
    # see Notes 3rd Nov. '21
    # -> corrected : see Notes 2th May '22
    sig_up  = group.err_up[i]
    sig_low = group.err_down[i]
    sig_sys = group.sys_error[i]
    sig_rnd = (sig_up+sig_low)/2.
    q = np.sqrt(sig_rnd**2 + sig_sys**2) - sig_rnd
    if prime=="up":
        return  sig_up  + q 
    elif prime=="low":
        return  sig_low + q
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
    warnings.warn("\nWARNING: we are here using the method B to combine groups. \
    See notes of 27thOct.\n",)
    series = copy.copy(series)
    #Set all the group in the series to have the same bin
    #comb_bins = get_bins(series,n_bins=150,return_bins=True) 
    #set all the series' group such that they are shifted and have their (shifted) density distribution
    # find the reference group (most accurate overall)
    ref_index = get_ref_index(series) # this, provided that each Group has its own tot_error,doens't change
    ref_G     = copy.copy(series[ref_index])
    ignore_G  = [ref_G.name]
    combined_groups=[ref_G]
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
        combined_groups.append(tns_ref)
        ref_G = combine_groups_methodB(ref_G,tns_ref)
        tension_series =[]
        for G in series:
            if G.name not in ignore_G and tau_G_methodB(G,ref_G)>=sigmathresh:
                tension_series.append(G)
    print("Combined:",ignore_G)
    ref_G.name+=r"\nCombined result with $\tau_{thresh}=$"+str(sigmathresh)
    ref_G.combined_names = ignore_G
    if return_combined_list:
        return ref_G,ignore_G,combined_groups
    else:
        return ref_G
