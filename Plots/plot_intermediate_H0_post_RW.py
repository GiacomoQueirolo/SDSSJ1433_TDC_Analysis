#!/usr/bin/env python
# coding: utf-8


"""
Plot to obtain intermediate H0 from single df of single filter model of lenstronomy and overplot the posterior
REWORKED
""";


import sys
import argparse
import matplotlib
import numpy as np
import pickle
import copy
import matplotlib.pyplot as plt

from Utils.tools import *
from TD_analysis import pycs_get_res
from Utils.get_res import *
from Utils.Dt_from_Df_reworked import *
from H0.H0_Combined_reworked import *
from Utils.statistical_tools import quantiles
from Plots.plotting_tools import base_colors

from H0.tools_H0 import *
from H0.H0_Combined_reworked import get_PH0




font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rc('figure',**{'figsize':(12,9)})
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)

if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################
    parser = argparse.ArgumentParser(description="Temporary plot of the posterior of H0 from the given lens models")
    parser.add_argument("-dtd","--dt_dataset",type=str, dest="dt_dataset", default="J1433_forcen",
                    help="Dataset name for the time delay (default: J1433_forcen)")
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                    help="Number of bins per dimension for the histog. sampling of the Df (be careful with it! too many bins can be catastrophic)")
    parser.add_argument("-dir","--directory",type=str, dest="dir", default=".",
                    help="In which directory to save it")
    parser.add_argument("-sl","--simple_legend",dest="simple_legend", default=False,action="store_true",
                        help="Draw a simplified legend with only the name of the filters")

    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args  = parser.parse_args()
    # bins for the fermat pot
    bins  = args.nbins 
    setting_names =  args.SETTING_FILES
    
    #Dt:
    dt_dataset=args.dt_dataset
    print(f"Default time delay result obtained from: {dt_dataset}")
    pycs_path   = "./my_pycs_scripts/"
    dt_comb_res = pycs_get_res.get_combined_res(dt_dataset,main_dir_path=pycs_path)
    kwargs_dt   = get_kwdt(dt_comb_res)
    backup_res_lnstr = "./backup_results/"
    savefig_path     = backup_res_lnstr+"/Post_H0/"+dir
    mkdir(savefig_path)
    
    H0s  = []
    PH0s = []
    f,ax = plt.subplots(1,1,figsize=(12,12))     
    for i,sets in enumerate(setting_names):
        print("setting:",sets)
        mcmc_fermat = get_mcmc_fermat(sets)
        mcmc_Df     = np.transpose(mcmc_fermat)[1:]-np.transpose(mcmc_fermat)[0]
        mcmc_c      = np.array(copy.deepcopy(mcmc_Df))
        mcmc_BC     = mcmc_c[1] - mcmc_c[0]  # BC = C - B = (C-A)-(B-A) = AC - AB
        mcmc_c[2]   = mcmc_BC
        mcmc_Df     = mcmc_c 
        print("WARNING: Given the low S/N of image D, I will discard here and instead consider Delta BC")    
        Post_Df,Post_Df_bins = np.histogramdd(mcmc_Df.T,bins=bins,density=True) 
        PH0,H0 = get_PH0(Df_dens=Post_Df,Df_dens_bins=Post_Df_bins,Dt_kw=kwargs_dt,setting=get_setting_module(sets,1))
        PH0s.append(PH0)
        H0s.append(H0)
        ax.scatter(H0,PH0,c=base_colors[i],label=strip_setting_name(sets,filter=simple_legend))
        h0_res,err_min,err_max= quantiles(PH0,H0,return_quantiles=False)
        yh0 = max(PH0)/2
        ax.errorbar(h0_res,yh0,yerr=None,xerr=[[err_min],[err_max]],fmt="r",capsize=4,c=base_colors[i])

    y_max = max([max(ph0) for ph0 in PH0s])
    ax.axvline(h0planck.H0,label="Planck",c="r",ls="--")
    ax.axvline(h0licow.H0,label="H0LiCOW",c="g",ls="--")
    ax.fill_between(np.linspace(h0planck-h0planck.sigma_min,h0planck+h0planck.sigma_max ) , -10, 10, color='red', alpha=0.2)
    ax.fill_between(np.linspace(h0licow-h0licow.sigma_min ,h0licow+h0licow.sigma_max) , -10, 10, color='green', alpha=0.2)
    plt.ylim(0,y_max*1.1)
    plt.xlabel("$H_0$[km/s/Mpc]")
    plt.ylabel("P($H_0$)")
    plt.legend()
    plt.title("Compare $H_0$ posterior from individual filter's lens models")
    f.savefig(savefig_path+"/compare_H0.png")
    print("Created "+savefig_path+"/compare_H0.png")
    
    res_H0 = [H0s,PH0s]
    with open(savefig_path+"/H0s.pkl","wb") as f:
        pickle.dump(res_H0,f)
    success(sys.argv[0])    
    

