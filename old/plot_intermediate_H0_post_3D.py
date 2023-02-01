#!/usr/bin/env python
# coding: utf-8

# In[80]:


"""
Plan to obtain intermediate H0 from single df of single filter model of lenstronomy and overplot the posterior
""";


# In[1]:


import os,sys
import argparse
import importlib
import matplotlib
import numpy as np
import corner,pickle
import json,copy,time
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from tools import *
import pycs_get_res
from get_res import *
from Dt_from_Df import *
from plotting_tools import base_colors


# In[ ]:


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rc('figure',**{'figsize':(12,9)})
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)


# In[ ]:


# Litt. Data
h0planck = 67.4 # km/s/Mpc
h0planck_err = 0.5
h0licow  = 73.3 # km/s/Mpc
h0licow_err = [1.8,1.7]


# In[ ]:


def get_bin_center(bins):
    if len(np.shape(bins))==1: #so it works for 1D too
        bins=[bins]
    bins_centers = []
    for i in range(len(bins)):
        bins_centers_i = []
        for j in range(len(bins[i])-1):
            bins_centers_i.append((bins[i][j]+bins[i][j+1])/2.)
        bins_centers.append(bins_centers_i)
    bins_centers =np.array(bins_centers)
    return bins_centers

def err_Df(err_dt,h0):
    # Convert sigma dt in sigma dphi 
    #-> simply the sqrt of the covariance
    cov_Df_ = cov_Df(err_dt**2,h0)
    err_Df_ = np.sqrt(cov_Df_)
    return err_Df_

def get_analytic_density_1D(mean,err,bins,norm=None):
    center_bins = get_bin_center(bins) 
    #grid_of_points = np.transpose(np.meshgrid(*center_bins))
    dens  = norm.pdf(center_bins,loc=mean,scale=cov)
    if type(norm) is float or type(norm) is int:
        dens/=norm
    elif type(norm) is bool:
        if norm:
            dens/=np.sum(dens)
    return dens


def get_PH0_1D(Dens_f,
            nd_bins,
            Dt_kw,
            H0=np.arange(35,100,.1)):#(Post_Df,Post_Df_bins,Dt_kw=dt_kws[dim_i])
    PH0 = []
    for h0 in H0:
        kwargs_df = {"mean":Df_XY(copy.deepcopy(kwargs_dt["mean"]),h0),
             "err":err_Df(copy.deepcopy(kwargs_dt["err"]),h0)}
        Dens_f_trsf = get_analytic_density_1D(bins=copy.deepcopy(nd_bins),**kwargs_df)
        Dens_tot    = copy.deepcopy(Dens_f)*Dens_f_trsf
        # the integration in this case is nothing else then the sum over every bin, ie:
        P_h0 = np.sum(Dens_tot)
        if np.isnan(P_h0):
            P_h0 = 0.
        PH0.append(P_h0)
    if np.sum(PH0)!=0:
        PH0 = PH0/np.sum(PH0)
    else:
        print("WARNING: sum of PH0 == 0")
    return np.array(PH0),H0


def quantiles_uncertainties(prob,sampling_prob,q=[0.16,.5,0.84],return_quantiles=False):
    qnt    = []
    integr = 0
    for i in range(len(prob)):
        integr+=prob[i]
        for qi in q:
            if integr-prob[i]<qi and integr>=qi:
                qnt.append(sampling_prob[i])
    if return_quantiles:
        return qnt
    else:
        res     = qnt[1]
        err_min = qnt[1]-qnt[0]
        err_max = qnt[2]-qnt[1]
        return (res,err_min,err_max)


# In[83]:


if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################
    parser = argparse.ArgumentParser(description="Temporary plot of the posterior of H0 from the given lens models")
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                    help="Number of bins per dimension for the histog. sampling of the Df (be careful with it! too many bins can be catastrophic)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    
    parser.add_argument('-ln','--lensname',help="Lensname for time delay dataset to consider",
                        dest="lensname", 
                        default="J1433",action="store")
    parser.add_argument('-dn','--dataname',action='store',default="forcen",
                        help="Dataname for time delay dataset to consider")

    args  = parser.parse_args()
    # bins for the fermat pot
    bins  = args.nbins 
    setting_names =  args.SETTING_FILES
    lensname = args.lensname
    dataname = args.dataname
    
    #Dt:
    print("Time delay result obtained from: "+lensname+"_"+dataname)
    pycs_path = "./my_pycs_scripts/"
    dt_comb_res = pycs_get_res.get_combined_res(lensname+"_"+dataname,main_dir_path=pycs_path)
    sys = dt_comb_res.error.sys
    rnd = dt_comb_res.error.rnd
    err_dt  = np.sqrt(sys**2 + rnd**2)
    dt_kws  = [{"mean":dt_comb_res.results[i],"err":err_dt[i]} for i in range(len(err_dt))]
    
    backup_res_lnstr = "./backup_results/"
    savefig_path     = backup_res_lnstr+"/Post_H0/"
    mkdir(savefig_path)
    
    mcmc_Dfs = []
    for i,sets in enumerate(setting_names):
        print("setting:",sets)
        mcmc_fermat = get_mcmc_fermat(sets)
        mcmc_Df     = np.transpose(mcmc_fermat)[1:]-np.transpose(mcmc_fermat)[0]
        mcmc_c    = np.array(copy.deepcopy(mcmc_Df))
        mcmc_BC   = mcmc_c[1] - mcmc_c[0]  # BC = C - B = (C-A)-(B-A) = AC - AB
        mcmc_c[2] = mcmc_BC
        mcmc_Df   = mcmc_c 
        if i==0:
            print("WARNING: Given the low S/N of image D, I will discard here and instead consider Delta BC")    
            lcs = ["AB","AC","BC"] 
        mcmc_Dfs.append(mcmc_Df)
        
    H0s  = []
    PH0s = []
    legend_elements = []
    f,axes = plt.subplots(2,2,figsize=(12,12)) 
    
    for dim_i in range(len(mcmc_Df)): # each dimension/quadrant
        if dim_i==0:
            ax = axes[0][0]
        elif dim_i==1:
            ax = axes[1][0]
        elif dim_i==2:
            ax = axes[1][1]
        else:
            raise RuntimeError("")
        y_max = 0
        for i,sets in enumerate(setting_names): #each setting
            Post_Df,Post_Df_bins = np.histogram(mcmc_Dfs[i][dim_i],bins=bins,density=True)     
            PH0,H0 = get_PH0_1D(Post_Df,Post_Df_bins,Dt_kw=dt_kws[dim_i])
            ax.scatter(H0,PH0,c=base_colors[i]) # label=strip_setting_name(sets)
            h0_res,err_min,err_max= quantiles_uncertainties(PH0,H0,return_quantiles=False)
            yh0 = max(PH0)/2
            ax.scatter(h0_res,yh0,c=base_colors[i])
            ax.errorbar(h0_res,yh0,yerr=None,xerr=[[err_min],[err_max]],capsize=4,fmt=base_colors[i])
            str_res = str(np.round(h0_res,2))+"$_{-"+str(np.round(err_min,2))+"}^{+"+str(np.round(err_max,2))+"}$"
            ax.text(h0_res-.5*len(str((np.round(h0_res,2)))),yh0*1.05,str_res,c=base_colors[i])
            if dim_i==0:
                legend_elements.append(Patch(facecolor=base_colors[i],label=strip_setting_name(sets)))
            if max(PH0)>y_max:
                y_max = max(PH0)
        ax.set_ylim(0,y_max*1.1)
        ax.axvline(h0planck,label="Planck",c="r",ls="--")
        ax.axvline(h0licow,label="H0LiCOW",c="g",ls="--")
        ax.fill_between(np.linspace(h0planck-h0planck_err,  h0planck+h0planck_err ), -10, 10, color='r', alpha=0.2)
        ax.fill_between(np.linspace(h0licow-h0licow_err[0], h0licow+h0licow_err[1]), -10, 10, color='g', alpha=0.2)
        ax.set_xlabel("$H_0$[km/s/Mpc] from "+lcs[dim_i])
        ax.set_ylabel("P($H_0$)")
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    plt.tight_layout()
    f.suo_title("Compare $H_0$ posterior from individual filter's lens models for different LCs combinations")
    f.savefig(savefig_path+"/compare_H0_3D.png")
    print("Created "+savefig_path+"/compare_H0_3D.png")
    
    res_H0 = [H0s,PH0s]
    with open(savefig_path+"/H0s.pkl","wb") as f:
        pickle.dump(res_H0,f)
    success(sys.argv[0])    

