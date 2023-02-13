#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# My version of params from kwargs to args 


# In[1]:


from lenstronomy.Sampling.parameters import Param
import importlib
import numpy as np
import sys
sys.path.append("..")
from Utils.tools import *

def my_kwarg2args(kwargs_res,setting_name,cut_mcmc=0):
    setting_module = get_setting_module(setting_name)
    setting        =     setting_module.setting()
    WS = check_if_WS(setting)
    
    if not check_if_CP(setting):
        lens_model_list = ['SIE']
    else:
        print("WARNING: Considering the PEMD profile for the main lens")
        lens_model_list = ['PEMD']
    lens_model_list= [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']
    
    # Lens light profile with perturber
    if setting.sub==False:
        light_model_list = ["SERSIC_ELLIPSE", "SERSIC","UNIFORM"]
    else:
        light_model_list = ["SERSIC","UNIFORM"]

    if hasattr(setting,"no_pert"):
        if setting.no_pert==True:
            light_model_list=["UNIFORM"]
    else:
        setting.no_pert = False
    # Source host gal
    if WS==False:
        source_model_list = ['SERSIC_ELLIPSE']
    # Source light model: point 
    point_source_list = ['LENSED_POSITION']

    kwargs_model = {'lens_model_list': lens_model_list,
                    'lens_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list,
                    'additional_images_list': [False], #list of bools (same length as point_source_type_list).
                    # If True, search for additional images of the same source is conducted.
    # 'fixed_magnification_list': [False],  # list of bools (same length as point_source_type_list).
                     'fixed_magnification_list': setting.fixed_mag,  # list of bools (same length as point_source_type_list).
                    #If True, magnification ratio of point sources is fixed to the one given by the lens model 
                    }
    if WS==False:
        kwargs_model['source_light_model_list'] = source_model_list
        
    if setting.sub==False:
        joint_lens_with_light=[[0,0,["center_x","center_y"]],[1,1,["center_x","center_y"]]]
    else:
        joint_lens_with_light=[[0,1,["center_x","center_y"]]]
        
    if setting.no_pert:
        joint_lens_with_light=[]
    
    kwargs_constraints = {'num_point_source_list': [4], #-> number of point source images
                               #'solver_type': 'PROFILE_SHEAR',  # 'PROFILE', \
                              'solver_type': 'NONE',
                              #'PROFILE_SHEAR', 'ELLIPSE', 'CENTER', 'NONE'
        #:param solver_type: string, option for specific solver type -> still not clear
                          'joint_lens_with_light':joint_lens_with_light
                         }
    if WS == False:
        kwargs_constraints['joint_source_with_point_source'] = [[0, 0]]
        
    kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens = setting.lens_params 
    kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps           = setting.ps_params
    kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light = setting.lens_light_params

    if WS==False:
        kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source = setting.source_params
        param_class = Param(kwargs_model, fixed_lens, fixed_source,fixed_lens_light, fixed_ps, [], 
                        [],kwargs_lower_lens, kwargs_lower_source,
                        kwargs_lower_lens_light, kwargs_lower_ps,[], 
                        [],kwargs_upper_lens, kwargs_upper_source, 
                        kwargs_upper_lens_light, kwargs_upper_ps,[],
                        [], **kwargs_constraints)
    else:
         
        param_class = Param(kwargs_model, fixed_lens, [],fixed_lens_light, fixed_ps, [], 
                        [],kwargs_lower_lens, [],kwargs_lower_lens_light, kwargs_lower_ps,[], 
                        [],kwargs_upper_lens, [], kwargs_upper_lens_light, kwargs_upper_ps,[],
                        [], **kwargs_constraints)
    #param_class.args2kwargs(bestfit_result, bijective=True)
    args_res = param_class.kwargs2args(**kwargs_res)
    return args_res


# In[ ]:


# Older version, same but probably less elastic
"""
import pickle
import numpy as np
import json

savemcmc_path="./"+backup_path+"/"+setting_name.replace("settings_","mcmc_")+"/"
savefig_path="./"+backup_path+"/"+setting_name.replace("settings_","")+"/" 
mcmc_file=savemcmc_path+setting_name.replace("settings","mcmc_smpl")+".json"
param_file=savemcmc_path+setting_name.replace("settings","mcmc_prm")+".dat"

cut_mcmc = 0

#MCMC sample
mcmc_file = open(mcmc_file, 'r')
samples_mcmc = np.array(json.load(mcmc_file))[cut_mcmc:]
mcmc_file.close()

#parameters' name
param_mcmc=[]
with open(param_file,"r") as param_file:
    param_mcmc=(param_file.readlines())
    
for i in range(len(param_mcmc)):
    param_mcmc[i]=param_mcmc[i].replace(",\n","")


with open(savefig_path+"/read_results.data","rb") as file:
    rres = pickle.load(file)
rres ,param_mcmc
rres_val =[]
count_ps_ra = 0
count_ps_dec = 0
for p in param_mcmc:
    prf=None
    for nm_prof in ["lens_light","source"]:
        if nm_prof in p:
            prf=nm_prof
    if prf==None:
        if "lens" in p:
            prf="lens"
        elif "image" in p:
            prf="ps"
    if prf!="ps":
        prm=None
        for nm_prm in ["theta_E","e1","e2","center_x","center_y","gamma_ext","psi_ext","R_sersic","n_sersic"]:
            if nm_prm in p:
                prm=nm_prm
        nm=None
        for n in range(4):
            if str(n) in p[-1]:
                nm=int(n)
        rres_val.append(rres["kwargs_"+prf][nm][prm])
    else:
        if "ra" in p:
            rres_val.append(rres["kwargs_"+prf][0]["ra_image"][count_ps_ra])  
            count_ps_ra+=1
        elif "dec" in p:
            rres_val.append(rres["kwargs_"+prf][0]["dec_image"][count_ps_dec])   
            count_ps_dec+=1


"""

