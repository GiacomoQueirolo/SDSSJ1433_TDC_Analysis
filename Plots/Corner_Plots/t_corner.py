#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Mainly for Dt_AD to check why the larger uncertainty and 
# to check if it actually depends on the center_y_lens0 parameter


# In[ ]:


import numpy as np
import os,sys
import corner
import matplotlib.pyplot as plt
from copy import copy
from astropy import units as u
import argparse

####################
from Utils.tools import *
from Utils.get_res import *
from Utils.Dt_from_Df_reworked import Dt_XY


# In[ ]:


param_mcmc_t = ["td_A","td_B","td_C","td_D"]
H0_planck = 67.4*u.km/u.second/u.megaparsec # km/s/Mpc

def t_corner(setting,cut_mcmc,save=False):
    setting = get_setting_module(setting).setting()
    savefig_path  = get_savefigpath(setting)
    #MCMC sample
    samples_mcmc = np.array(get_mcmc_smpl(setting,backup_path)[cut_mcmc:])
    #MCMC sample for Fermat pot
    #NOTE: they are ordered A,B,C,D
    mcmc_fermat = np.array(get_mcmc_fermat(setting,backup_path)[cut_mcmc:])
    #parameters' name
    param_mcmc =get_mcmc_prm(setting,backup_path)
    # Time delay behaviour
    ########################
    mcmc_t = Dt_XY(Df_XY=mcmc_fermat,H0=H0_planck,z_l=setting.z_lens,z_s=setting.z_source)
    # order A, B, C, D
    new_samples = copy(samples_mcmc.T)
    for t in mcmc_t.T:
        new_samples = np.array([*new_samples,t])
    param_mcmc = param_mcmc + param_mcmc_t
    red_samples = []
    red_param = []
    for sample_i,param_i in zip(new_samples,param_mcmc):
        if "td_" in param_i or "center_" in param_i or "ra_" in param_i or "dec_" in param_i:
            red_samples.append(sample_i)
            prm_i = param_i
            if "center_" in param_i:
                prm_i = setting.str_param(param_i)
            red_param.append(prm_i)
    red_samples = np.array(red_samples).T
    
    mcmc_t  = np.array(mcmc_t)
    Delta_T = np.transpose(mcmc_t.T[1:]-mcmc_t.T[0]).tolist()
    labels_DT  = ["$\Delta t_{A-B}$","$\Delta t_{A-C}$", "$\Delta t_{A-D}$"]
    if save:
        plot= corner.corner(red_samples, labels=red_param,show_titles=True)
        plot.savefig(savefig_path+"/t_corner_red.png")

        # consider also the Dt plot (alone, but now with real Deltas)
        plot    = corner.corner(np.array(Delta_T), labels=labels_DT,show_titles=True)
        #plot.title("$H_0^{Planck} = $"+str(H0_planck))
        plot.savefig(savefig_path+"/Dt_corner.png")
    else:
        return red_samples,red_param,Delta_T,labels_DT


# In[ ]:


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Produces the MCMC posterior and corner plots of SOME parameters, Dt in particular")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument('SETTING_FILES',nargs="+",default=[], help="setting file(s) to consider")
    args = parser.parse_args()
    settings = [get_setting_module(sett).setting() for sett in args.SETTING_FILES]
    cut_mcmc = args.cut_mcmc
    ##########################
    present_program(sys.argv[0])
    ##########################
    backup_path="./backup_results/"
    for setts in settings:
        if len(settings)>1:
            print("Setting "+get_setting_name(setts)) 
        t_corner(setts,cut_mcmc,save=True)
    success(sys.argv[0])

