#!/usr/bin/env python
# coding: utf-8

# # Compute the fermat potential prior

# copy from Df_pot_prior, but more similar to prior_fermat_pot.ipynb
# basically instead of creating it from the setting file, it is based on the mcmc results and then we create from that the uniform distribution for the prior of the lens parameters

# In[2]:


# import of standard python libraries
import copy
import json
import corner
import os,sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from lenstronomy.LensModel.lens_model import LensModel

#%matplotlib inline
from Utils.tools import *
from Posterior_analysis.mag_remastered import get_mag_mcmc
from Posterior_analysis.fermat_pot_analysis import gen_mcmc_Df
from Utils.order_images import get_new_image_order
from Data.Param import get_prm_list
# In[1]:


# This parts give the mcmc of the fermat potential given the samples_mcmc prior


def get_mcmc_Df_prior(mcmc_prior,setting,Df_boundaries=None,threshold_mcmc_points=1000,previous_mcmc=None,Ntot_previous=0): 
    # mcmc_prior shape = param, n_points
    # setting = setting module
    # Df_boundaries = [[min_AB,max_AB],[min_A..],..]
    # threshold_mcmc_points = int, minumum number of poiints in the boundaries
    # previous_mcmc = previous points to append to the sampling
    mcmc_prior_Df = gen_mcmc_Df(mcmc_prior,setting)  #D_AB, D_AC, D_AD, meaning Df_i - Df_A
    if previous_mcmc is not None:
        mcmc_prior_Df=np.vstack((mcmc_prior_Df,previous_mcmc))    
    # cut points out of boundaries
    success = [True, "SUCCESS"]
    Ntot = Ntot_previous + len(mcmc_prior_Df)
    if Df_boundaries is not None:
        excluded_points=0
        mcmc_prior_Df_cut = []
        for mcmcdfi in mcmc_prior_Df:
            if not is_in_boundaries(mcmcdfi,Df_boundaries):
                excluded_points+=1
            else:
                mcmc_prior_Df_cut.append(mcmcdfi)
        # double check
        # for some reason necessary to eliminate 1 single point
        # that appeared where it was not supposed to
        for i in range(len(mcmc_prior_Df_cut)):
            if not is_in_boundaries(mcmc_prior_Df_cut[i],Df_boundaries):
                del mcmc_prior_Df_cut[i]
        if len(mcmc_prior_Df_cut)<threshold_mcmc_points:
            success_string = "The excluded points are {}, leaving only {} points in the final fermat MCMC. Run with larger n_points".format(excluded_points,len(mcmc_prior_Df_cut))
            success=[False,success_string]
        return mcmc_prior_Df_cut,success,Ntot
        
        #return mcmc_prior_Df_cut
    
    return mcmc_prior_Df,success,Ntot # shape: Df_ai, mcmc_steps
    #plot = corner.corner(mcmc_Df, labels=labels_Df, show_titles=True)
    #plot.savefig(savefig_path+"test_uniform_df.png")
    #print("test result in "+savefig_path+"test_uniform_df.png")


# In[ ]:


def get_mcmc_mag_prior(mcmc_prior,setting,mag_boundaries=None,threshold_mcmc_points=1000,previous_mcmc=None,Ntot_previous=0): 
    # mcmc_prior shape = param, n_points
    # setting = setting module
    # mag_boundaries = [[min_AB,max_AB],[min_A..],..]
    # threshold_mcmc_points = int, minumum number of poiints in the boundaries
    # previous_mcmc = previous points to append to the sampling
    CP          = check_if_CP(setting)
    param_mcmc  = get_prm_list(setting)
    mcmc_mag    = get_mag_mcmc(mcmc_prior,param_mcmc,setting,CP)
    #I want to obtain the correct image order
    #########################################
    new_order = get_new_image_order(setting,mcmc_prior,starting_from_A=False)
    tmp_mcmc  = np.transpose(mcmc_mag) 
    mcmc_prior_mag = np.transpose([tmp_mcmc[i] for i in new_order]).tolist() #mu_AB, mu_AC, mu_AD, meaning mu_i / mu_A
    
    if previous_mcmc is not None:
        mcmc_prior_mag=np.vstack((mcmc_prior_mag,previous_mcmc))    
    # cut points out of boundaries
    success = [True, "SUCCESS"]
    Ntot = Ntot_previous + len(mcmc_prior_mag)
    if mag_boundaries is not None:
        excluded_points=0
        mcmc_prior_mag_cut = []
        for mcmc_magi in mcmc_prior_mag:
            if not is_in_boundaries(mcmc_magi,mag_boundaries):
                excluded_points+=1
            else:
                mcmc_prior_mag_cut.append(mcmc_magi)
        # double check
        # for some reason necessary to eliminate 1 single point
        # that appeared where it was not supposed to
        for i in range(len(mcmc_prior_mag_cut)):
            if not is_in_boundaries(mcmc_prior_mag_cut[i],mag_boundaries):
                del mcmc_prior_mag_cut[i]
        if len(mcmc_prior_mag_cut)<threshold_mcmc_points:
            success_string = "The excluded points are {}, leaving only {} points in the final mag MCMC. Run with larger n_points".format(excluded_points,len(mcmc_prior_mag_cut))
            success=[False,success_string]
        return mcmc_prior_mag_cut,success,Ntot
    
    return mcmc_prior_mag,success,Ntot # shape: mag_ai, mcmc_steps


# In[1]:


"""# This parts give the mcmc of the mag ratio the samples_mcmc prior
def get_mcmc_mag_prior(mcmc_prior,param_mcmc,setting,mag_boundaries=None,threshold_mcmc_points=100): 
    # mcmc_prior shape = param, n_points
    # param_mcmc = list of name of params
    # setting = setting module
    # mag_boundaries = [[min_AB,max_AB],[min_A..],..]
    
    CP= check_if_CP(setting)
    mcmc_mag = get_mag_mcmc(mcmc_prior,param_mcmc,setting,CP)
  
    #I want to obtain the correct image order
    #########################################
    new_order = get_new_image_order(setting,mcmc_prior,starting_from_A=False)
    tmp_mcmc = np.transpose(mcmc_mag) 
    mcmc_prior_mag_rt = np.transpose([tmp_mcmc[i] for i in new_order]).tolist() #mu_AB, mu_AC, mu_AD, meaning mu_i / mu_A
    
    # cut points out of boundaries
    if mag_boundaries is not None:
        excluded_points=0
        mcmc_prior_mag_cut = []
        for mcmcMui in mcmc_prior_mag_rt:
            if not is_in_boundaries(mcmcMui,mag_boundaries):
                excluded_points+=1
            else:
                mcmc_prior_mag_cut.append(mcmcMui)
        if len(mcmc_prior_mag_cut)<threshold_mcmc_points:
            raise UserWarning("The excluded points are",excluded_points,", leaving only",len(mcmc_prior_mag_cut)," points in the final Mag MCMC. Run with larger n_points")
        
        # double check
        # for some reason necessary to eliminate 1 single point
        # that appeared where it was not supposed to
        for i in range(len(mcmc_prior_mag_cut)):
            if not is_in_boundaries(mcmc_prior_mag_cut[i],mag_boundaries):
                del mcmc_prior_mag_cut[i]
        return mcmc_prior_mag_cut
    
    return mcmc_prior_mag_rt # shape: mag_ai, mcmc_steps
"""


# In[ ]:


def is_in_boundaries(dfi,Df_boundaries):
    for IJ,Dfi_IJ in enumerate(dfi):
        if not Dfi_IJ>=Df_boundaries[IJ][0] or not Dfi_IJ<Df_boundaries[IJ][1]:            
            return False
    return True


# In[ ]:


def get_prior(setting,npoints):
    param_mcmc = get_prm_list(setting)
    n_images   = count_images(param_mcmc)
    kw_models  = get_kwargs_to_model(setting)
    mcmc_prior = []
    for ip,param in enumerate(param_mcmc):
        if "image" in param:
            par_min,par_max = get_boundary_param(param,kw_models,ip,len(param_mcmc),n_images)
        else:
            par_min,par_max = get_boundary_param(param,kw_models)
        mcmc_prior.append(np.random.uniform(par_min,par_max,npoints))
    mcmc_prior = np.transpose(mcmc_prior).tolist()
    return mcmc_prior
    

def count_images(params):
    n_im_ra  = params.count("ra_image")
    n_im_dec = params.count("dec_image")
    """
    for p in params:
        if "ra_image" in p:
            n_im_ra+=1
        elif "dec_image" in p:
            n_im_dec+=1
    """
    if n_im_dec!=n_im_ra:
        raise RuntimeError("Number of coordinates for images in parameters not matching")
    return n_im_ra

def get_kwargs_to_model(setting):
    setting = get_setting_module(setting).setting()
    kw_to_model= {}
    #those are always present
    kw_to_model["lens_params"]       = setting.lens_params
    kw_to_model["ps_params"]         = setting.ps_params
    kw_to_model["lens_light_params"] = setting.lens_light_params  # at least the uniform
    try:
        kw_to_model["source_params"] = setting.source_params
    except:
        pass
    return kw_to_model

def get_boundary_param(param,kw_models,ip=None,nparam=None,n_images=4):
    kw_model,split_name = kw_model_from_param(kw_models,param)
    lower_kw = kw_model[3]
    upper_kw = kw_model[4]
    if split_name=="_image" and ip is None:
        raise RuntimeError("Give me the index of the parameter for the ps image boundaries")
    elif split_name=="_image" and nparam is None:
        raise RuntimeError("Give me the lenght of the parameters for the ps image boundaries")
    elif split_name=="_image":
        nimg = index_image(param,ip,nparam,n_images)
        par_min,par_max = lower_kw[0][param][nimg],upper_kw[0][param][nimg]
    else:
        index_model = int(param[-1])
        key_prm = param.split(split_name)[0] 
        par_min = lower_kw[index_model][key_prm]
        par_max = upper_kw[index_model][key_prm]
    return par_min,par_max

def kw_model_from_param(kw_models,param):
    for word in ["theta_E","e1","e2","gamma","gamma_ext","psi_ext"]:
        if word+"_lens" in param:
            return kw_models["lens_params"],"_lens"
    for word in ["ra","dec"]:
        if word+"_image" in param:
            return kw_models["ps_params"],"_image"
    if "lens_light" in param: 
        # this should catch if there are "center_x/y_lens_light"
        # so that i consider the correct prior
        return kw_models["lens_light_params"],"_lens_light"
    if "source_light" in param:
        try:
            return kw_models["source_params"],"_source_light"
        except KeyError:
            raise KeyError("The given setting file do not have source params, but the resulting mcmc have them")
    
    if "center_x_lens" in param or "center_y_lens" in param:
        return kw_models["lens_params"],"_lens"
    # if it doesn't manage to find anything, creates an error
    return None
    
def index_image(param,ip,nparam,n_images=4):
    # get the index of the image
    n_starting = nparam - 2*n_images #n_images*ra,n_images*dec
    if "ra_" in param:
        return int(ip-n_starting)
    elif "dec_" in param:
        return int(ip-n_starting-n_images)
    else:
        raise RuntimeError("These params should be image")


# In[ ]:


def Df_prior(setting,npoints=10000,Df_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_Df_prior"):
    setting_name = get_setting_name(setting)
    setting = get_setting_module(setting_name).setting()
    mcmc_Df_prior_ = None 
    while True:
        mcmc_prior  = get_prior(setting,npoints)
        mcmc_Df_prior,success,Ntot  = get_mcmc_Df_prior(mcmc_prior,setting,Df_boundaries,previous_mcmc=mcmc_Df_prior_)
        if success[0]:
            break
        else:
            if np.shape(mcmc_Df_prior)!=np.shape([]):
                mcmc_Df_prior_     = copy.copy(mcmc_Df_prior)
            print(success[1])
            print("Increasing n_points of 70000")
            npoints += 70000
    print("Df prior Done")
    if save_mcmc:
        mcmc_path = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
        mcmc_prior_path = mcmc_path+"/"+output_name+".json"
        with open(mcmc_prior_path,"w") as f:
            json.dump(np.array(mcmc_Df_prior).tolist(),f)
    return mcmc_Df_prior,Ntot


# In[ ]:


def Df_prior_ABC(setting,npoints=10000,Df_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_Df_prior_ABC"):
    setting_name = get_setting_name(setting)
    # Consider ABC but not D
    mcmc_Df_prior,Ntot = Df_prior(setting_name,npoints,Df_boundaries=Df_boundaries,save_mcmc=False)
    
    print("WARNING: Only considering ABC")    
    mcmc_c    = np.array(copy.deepcopy(mcmc_Df_prior))
    # BC = C - B = (C-A)-(B-A) = AC - AB
    mcmc_BC   = mcmc_c[1]- mcmc_c[0]
    mcmc_c[2] = mcmc_BC
    mcmc_Df_prior_ABC = mcmc_c
    
    if save_mcmc:
        mcmc_path = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
        mcmc_prior_path = mcmc_path+"/"+output_name+".json"
        with open(mcmc_prior_path,"w") as f:
            json.dump(mcmc_Df_prior_ABC.tolist(),f)
            
    return mcmc_Df_prior_ABC,Ntot


# In[ ]:


def mag_prior(setting,npoints=10000,mag_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_mag_rt_prior"):
    setting = get_setting_module(setting,1)
    mcmc_mag_prior_ = None 
    while True:
        mcmc_prior = get_prior(setting,npoints)
        mcmc_mag_prior,success,Ntot  = get_mcmc_mag_prior(mcmc_prior,setting,mag_boundaries,previous_mcmc=mcmc_mag_prior_)
        if success[0]:
            break
        else:
            if np.shape(mcmc_mag_prior)!=np.shape([]):
                mcmc_mag_prior_     = copy.copy(mcmc_mag_prior)
            print(success[1])
            print("Increasing n_points of 70000")
            npoints += 70000
    print("Mag prior Done")
    if save_mcmc:
        mcmc_path       = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
        mcmc_prior_path = mcmc_path+"/"+output_name+".json"
        with open(mcmc_prior_path,"w") as f:
            json.dump(np.array(mcmc_mag_prior).tolist(),f)
    return mcmc_mag_prior,Ntot


def mag_prior_ABC(setting,npoints=10000,mag_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_mag_rt_prior_ABC"):
    # Consider ABC but not D
    setting_name = get_setting_name(setting)
    mcmc_mag_prior,Ntot = mag_prior(setting_name,npoints,mag_boundaries=mag_boundaries,save_mcmc=False)
    print("WARNING: Only considering ABC")    
    mcmc_c    = np.array(copy.deepcopy(mcmc_mag_prior))
    # BC = C / B = (C/A)/(B/A) = AC/AB
    mcmc_BC   = mcmc_c[1]/mcmc_c[0]
    mcmc_c[2] = mcmc_BC
    mcmc_mag_prior_ABC = mcmc_c
    
    if save_mcmc:
        mcmc_path = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
        mcmc_prior_path = mcmc_path+"/"+output_name+".json"
        with open(mcmc_prior_path,"w") as f:
            json.dump(mcmc_mag_prior_ABC.tolist(),f)
            
    return mcmc_mag_prior_ABC,Ntot


# In[ ]:


"""
def mag_prior(setting,npoints=10000,mag_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_mag_rt_prior"):
    setting_name = get_setting_name(setting)
    mcmc_path = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
    mcmc_prior_path = mcmc_path+"/"+output_name+".json"

    setting = get_setting_module(setting_name).setting()
    
    while True:
        mcmc_prior,param_mcmc = get_prior(setting,npoints)
        try:
            mcmc_mag_prior  = get_mcmc_mag_prior(mcmc_prior,param_mcmc,setting,mag_boundaries)
            break
        except UserWarning:
            print("Increasing n_points of 70000")
            npoints       += 70000
    print("Mag prior Done")
    if save_mcmc:
        with open(mcmc_prior_path,"w") as f:
            json.dump(np.array(mcmc_mag_prior).tolist(),f)
    return mcmc_mag_prior

def mag_prior_ABC(setting,npoints=10000,mag_boundaries=None,save_mcmc=False,backup_path = "./backup_results/",output_name="mcmc_mag_rt_prior_ABC"):
    # Consider ABC but not D
    setting_name = get_setting_name(setting)
    mcmc_path = backup_path+"/"+setting_name.replace("settings","mcmc").replace(".py","")
    mcmc_prior_path = mcmc_path+"/"+output_name+".json"
    mcmc_mag_prior = mag_prior(setting_name,npoints,mag_boundaries=mag_boundaries,save_mcmc=False,backup_path=backup_path,output_name=output_name)
    
    print("WARNING: Only considering ABC")    
    mcmc_c    = np.array(copy.deepcopy(mcmc_mag_prior))
    # BC = C / B = (C/A)/(B/A) = AC/AB
    mcmc_BC   = mcmc_c[1]/mcmc_c[0]
    mcmc_c[2] = mcmc_BC
    mcmc_mag_prior_ABC = mcmc_c
    
    if save_mcmc:
        with open(mcmc_prior_path,"w") as f:
            json.dump(mcmc_mag_prior_ABC.tolist(),f)
    return mcmc_mag_prior_ABC

"""


# In[ ]:



