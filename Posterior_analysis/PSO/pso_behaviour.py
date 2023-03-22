#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#With argparse which should help make it easier
.ipynb_checkpoints/Notes-checkpoint.ipynb

# In[1]:

import os,sys
import corner
import argparse
import importlib
import numpy as np
from scipy import stats
from datetime import datetime
import matplotlib.pyplot as plt


from Utils.tools import *
from Utils.get_res import *
from Plots.plotting_tools import base_colors

def pso_behaviour(ax, pso_pos, param_name, col="b", num_average=1):
    """
    plots the pso behaviour and looks for convergence of the chain
    :param samples_mcmc: parameters sampled 2d numpy array
    :param param_mcmc: list of parameters
    :param num_average: number of samples to average (should coincide with the number of samples in the emcee process)
    :return:
    """
    num_samples = len(pso_pos)
    num_average = int(num_average)
    n_points    = int((num_samples - num_samples % num_average) / num_average)

    samples = pso_pos
    sliced_sample = samples[:int(n_points * num_average)].reshape(n_points, num_average)
    samples_averaged = np.average(sliced_sample, axis=1)
    end_point = np.mean(samples_averaged)
    samples_renormed = (samples_averaged - end_point) / np.std(samples_averaged)
    x=np.arange(0,len(samples_renormed))
    ax.plot(x,samples_renormed, label=param_name,color=col)
    mn = min(samples_renormed)
    fac_mn = 1-0.002
    if mn<0:
        fac_mn = 1+0.002
    mx = max(samples_renormed)
    fac_mx = 1+0.002
    if mx<0:
        fac_mx = 1-0.002
    ax.set_ylim(mn*fac_mn,mx*fac_mx)
    ax.set_xlabel("PSO Steps")
    ax.set_ylabel("Norm. PSO position")
    ax.legend()
    return ax

"""def pso_simil_behaviour(ax, pso_pos, param_name,col="b", num_average=1,boundary=None):
"""
    #plots the pso behaviour and looks for convergence of the chain
    #:param pso_pos: parameters positions
    #:param param_mcmc: list of parameters
    #:param num_average: number of samples to average 
    #    (should coincide with the number of samples in the emcee process)
    #:return:
"""
    num_samples = len(pso_pos)
    num_average = int(num_average)
    n_points    = int((num_samples - num_samples % num_average) / num_average)

    samples = pso_pos
    sliced_sample    = samples[:int(n_points * num_average)].reshape(n_points, num_average)
    samples_averaged = np.average(sliced_sample, axis=1)
    
    x=np.arange(0,len(samples_averaged))
        
    
    ax.plot(x,samples_averaged,label=param_name,color=col)
    mn = min(samples_averaged)
    fac_mn = 1-0.002
    if mn<0:
        fac_mn = 1+0.002
    mx = max(samples_averaged)
    fac_mx = 1+0.002
    if mx<0:
        fac_mx = 1-0.002
    ax.set_ylim(mn*fac_mn,mx*fac_mx)
    
    #MOD_BOUNDARIES
    if boundary != None and "ra_" not in param_name and "dec_" not in param_name: 
        lim_up   = max(samples_averaged)
        lim_dwn  = min(samples_averaged)
        diff_lim = lim_up-lim_dwn
        lim_up  +=diff_lim/14.
        lim_dwn -=diff_lim/14.
        if min(boundary)>=(lim_dwn-4*diff_lim/14.):
            ax.hlines(boundary[0],0,len(samples_averaged),colors="r",linestyles='dashed',label="lower bound")
        elif max(boundary)<=(lim_up+4*diff_lim/14.):
            ax.hlines(boundary[1],0,len(samples_averaged),colors="r",linestyles='dashed',label="upper bound")
    ax.legend()
    ax.set_xlabel("PSO Steps")
    ax.set_ylabel("PSO position")
    return ax"""

def pso_sim_behaviour(ax, pso_pos, param_name,col="b", boundary=None):
    """
    plots the pso behaviour and looks for convergence of the chain
    :param pso_pos: parameters positions
    :param param_mcmc: list of parameters
    :param num_average: number of samples to average 
        (should coincide with the number of samples in the emcee process)
    :return:
    """
    
    x=np.arange(0,len(pso_pos))

    ax.plot(x,pso_pos,label=param_name,color=col)
    mn = min(pso_pos)
    fac_mn = 1-0.002
    if mn<0:
        fac_mn = 1+0.002
    mx = max(pso_pos)
    fac_mx = 1+0.002
    if mx<0:
        fac_mx = 1-0.002
    ax.set_ylim(mn*fac_mn,mx*fac_mx)
    
    #MOD_BOUNDARIES
    if boundary != None and "ra_" not in param_name and "dec_" not in param_name: 
        lim_up   = max(pso_pos)
        lim_dwn  = min(pso_pos)
        diff_lim = lim_up-lim_dwn
        lim_up  +=diff_lim/14.
        lim_dwn -=diff_lim/14.
        if min(boundary)>=(lim_dwn-4*diff_lim/14.):
            ax.hlines(boundary[0],0,len(pso_pos),colors="r",linestyles='dashed',label="lower bound")
        elif max(boundary)<=(lim_up+4*diff_lim/14.):
            ax.hlines(boundary[1],0,len(pso_pos),colors="r",linestyles='dashed',label="upper bound")
    ax.legend()
    
    ax.set_xlabel("PSO Steps")
    ax.set_ylabel("PSO position")
    return ax


def find_bound(param_name, setting,num=None):
    if "lens_light" in param_name:
        lens_light_params = setting.lens_light_params 
        n_lens = int(param_name[-1])
        lower  = lens_light_params[3][n_lens]       
        upper  = lens_light_params[4][n_lens]
        prm    = param_name.replace("_lens_light"+str(n_lens),"")
    elif "lens" in param_name:
        lens_params = setting.lens_params 
        n_lens = int(param_name[-1])
        lower  = lens_params[3][n_lens]       
        upper  = lens_params[4][n_lens]
        prm    = param_name.replace("_lens"+str(n_lens),"")
    elif "source_light0" in param_name:
        source_params = setting.source_params
        lower = source_params[3][0]       
        upper = source_params[4][0]        
        prm   = param_name.replace("_source_light0","")
    elif "image" in param_name:
        if "ra_" in param_name:
            radec="ra_image"
        else:
            radec="dec_image"
        ps_params = setting.ps_params
        lower     = ps_params[3][0][radec]
        upper     = ps_params[4][0][radec]
        prm       = num
    else:
        raise ValueError("Unrecognised parameter name:"+param_name)
    return lower[prm],upper[prm]


# In[ ]:


def plot_pso_chain(pso_chain,setting,savefig_path=None):
    _, steps, prm =  pso_chain
    lkl, pos, vel = steps
    pos_pso = np.transpose(pos)
    j=0
    n_ra=0
    n_dec=0
    axes=[]
    for i in range(len(prm)):  
        f, ax = plt.subplots(1, 1, figsize=(18, 6))
        param_i  = prm[i] 
        sample_t = pos_pso[i]
        if param_i=="dec_image" or param_i=="ra_image":
            param_i += "_"+str(i)
            udm      = "[\"]"
            if param_i=="dec_image":
                n_dec+=1
            else:
                n_ra+=1
        else:
            param_i,udm = setting.str_param(prm[i])
        if j>= len(base_colors):
            j=0
        col=base_colors[j]

        #MOD_BOUNDARIES
        if "ra_image" in param_i:
            min_bound,max_bound = find_bound(param_mcmc[i], setting,num=n_ra-1)
        elif "dec_image" in param_i:
            min_bound,max_bound = find_bound(param_mcmc[i], setting,num=n_dec-1)
        else:
            min_bound,max_bound = find_bound(param_mcmc[i], setting)


        f, ax = plt.subplots(1, 1, figsize=(18, 6))
        pso_sim_behaviour(ax=ax,pso_pos=sample_t,param_name=param_i, 
                                     boundary=[min_bound,max_bound],col=col)

        del_ax = ax.plot([],[])
        del_ax[0].set_color("r")
        title_mcmc_behaviour = f'PSO behaviour for {param_i} {udm}'
        ax.set_title(title_mcmc_behaviour)
        if savefig_path:
            plt.savefig(f'{savefig_path}/PSO_behaviour_{param_mcmc[i]}.png')
            plt.close()
        axes.append(ax)
        j+=1
    return axes


# In[ ]:


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the PSO behaviour plots of ALL parameter separately for the given filter model")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

    args     = parser.parse_args()
    present_program(sys.argv[0])
    setting_names = [st.replace(".py","") for st in args.SETTING_FILES]
    backup_path   = "backup_results/" 
    for setting_name in setting_names:
        if len(setting_names)>1:
            print_setting(setting_name)
        setting       = get_setting_module(setting_name).setting()
        savefig_path  = get_savefigpath(setting_name,backup_path)+"/PSO_bhv/"
        mkdir(savefig_path)
        savemcmc_path = get_savemcmcpath(setting_name,backup_path) 
        
        param_mcmc    = get_mcmc_prm(setting_name,backup_path)
        pso_file      = save_json_name(setting_name,path=savemcmc_path,filename="pso")
        pso_chain     = load_whatever(pso_file)
        plot_pso_chain(pso_chain=pso_chain,setting=setting,savefig_path=savefig_path)
    success(sys.argv[0])

