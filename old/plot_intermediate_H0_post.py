#!/usr/bin/env python
# coding: utf-8

# In[80]:


"""
Plan to obtain intermediate H0 from single df of single filter model of lenstronomy and overplot the posterior
""";


# In[81]:


import os,sys
import argparse
import matplotlib
import numpy as np
import corner,pickle
import json,copy,time
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

from tools import *
import pycs_get_res
from get_res import *
from Dt_from_Df import *
from plotting_tools import base_colors


# In[82]:


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)
matplotlib.rc('figure',**{'figsize':(12,9)})
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)


# In[92]:


def get_bin_center(bins):
    bins_centers = []
    for i in range(len(bins)):
        bins_centers_i = []
        for j in range(len(bins[i])-1):
            bins_centers_i.append((bins[i][j]+bins[i][j+1])/2.)
        bins_centers.append(bins_centers_i)
    bins_centers =np.array(bins_centers)
    return bins_centers

def get_analytic_density(mean,cov,bins,norm=1):
    center_bins = get_bin_center(bins) 
    grid_of_points = np.transpose(np.meshgrid(*center_bins))
    dens  = multivariate_normal.pdf(grid_of_points,mean,cov)*norm
    return dens


def get_normalisation(in_bins,Dt_kw,H0_range):
    bins = np.array(copy.deepcopy(in_bins))
    dH0 = H0_range[1]-H0_range[0] #is this necessary?
    kwargs_df = {"mean":Df_XY(copy.deepcopy(Dt_kw["mean"]),H0_range[0]),
     "cov":cov_Df(copy.deepcopy(Dt_kw["cov"]),H0_range[0])}
    norm = get_analytic_density(bins=bins,**kwargs_df)
    for i in range(1,len(H0_range)):
        kwargs_df = {"mean":Df_XY(copy.deepcopy(Dt_kw["mean"]),H0_range[i]),
             "cov":cov_Df(copy.deepcopy(Dt_kw["cov"]),H0_range[i])}
        norm+= get_analytic_density(bins=bins,**kwargs_df)*dH0
    print("normalisation done")
    return norm
def get_bin_vol(in_bins):
    bins =copy.deepcopy(in_bins)
    dim = np.shape(bins)[0]
    print(dim)
    bin_vol = 1
    for d in range(dim):
        bin_lenght=bins[d][1]-bins[d][0] 
        for n in range(5):
            rnd_index = np.random.randint(2,len(bins[d]))
            test_lenght =bins[d][rnd_index]-bins[d][rnd_index-1]
            if np.abs(bin_lenght-test_lenght)>1e-7:
                raise RuntimeError("Not all bins have the same size")
        bin_vol *= bin_lenght
    # we get a "volume" (vol in 3D, area in 2D, lenght in 1D) 
    # assuming that all bins are the same
    return bin_vol


# In[91]:


def get_PH0(Dens_f,
            nd_bins,
            Dt_kw,
            H0=np.arange(50,100,.1)):    
    PH0 = []
    for h0 in H0:
        kwargs_df = {"mean":Df_XY(copy.deepcopy(kwargs_dt["mean"]),h0),
             "cov":cov_Df(copy.deepcopy(kwargs_dt["cov"]),h0)}
        Dens_f_trsf = get_analytic_density(bins=copy.deepcopy(nd_bins),**kwargs_df)
        Dens_tot = copy.deepcopy(Dens_f)*Dens_f_trsf
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

    
# Litt. Data
h0planck = 67.4 # km/s/Mpc
h0planck_err = 0.5
h0licow  = 73.3 # km/s/Mpc
h0licow_err = [1.8,1.7]


# In[83]:


if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################
    parser = argparse.ArgumentParser(description="Temporary plot of the posterior of H0 from the given lens models")
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                    help="Number of bins per dimension for the histog. sampling of the Df (be careful with it! too many bins can be catastrophic)")
    parser.add_argument('-ln','--lensname',help="Lensname for time delay dataset to consider",
                        dest="lensname", 
                        default="J1433",action="store")
    parser.add_argument('-dn','--dataname',action='store',default="forcen",
                        help="Dataname for time delay dataset to consider")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
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
    cov_sys =  (np.array(dt_comb_res.error.sys)**2)*np.identity(len(dt_comb_res.error.sys))
    cov_rnd = np.cov(dt_comb_res.err_distr)
    cov_dt  = cov_sys + cov_rnd
    kwargs_dt = {"mean":dt_comb_res.results,"cov":cov_dt}
    
    backup_res_lnstr = "./backup_results/"
    savefig_path     = backup_res_lnstr+"/Post_H0/"
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
        
        PH0,H0 = get_PH0(Post_Df,Post_Df_bins,Dt_kw=kwargs_dt)
        PH0s.append(PH0)
        H0s.append(H0)
        ax.scatter(H0,PH0,c=base_colors[i],label=strip_setting_name(sets))
        h0_res,err_min,err_max= quantiles_uncertainties(PH0,H0,return_quantiles=False)
        yh0 = max(PH0)/2
        str_res = str(np.round(h0_res,2))+"$_{-"+str(np.round(err_min,2))+"}^{+"+str(np.round(err_max,2))+"}$"
        ax.errorbar(h0_res,yh0,yerr=None,xerr=[[err_min],[err_max]],fmt=base_colors[i],capsize=4)
        ax.scatter(h0_res,yh0,c=base_colors[i],fmt=".")
        ax.text(h0_res-.5*len(str((np.round(h0_res,2)))),yh0*1.05,str_res,c=base_colors[i])


    y_max = max([max(ph0) for ph0 in PH0s])
    ax.axvline(h0planck,label="Planck",c="r",ls="--")
    ax.axvline(h0licow,label="H0LiCOW",c="g",ls="--")
    ax.fill_between(np.linspace(h0planck-h0planck_err ,h0planck+h0planck_err ) , -10, 10, color='red', alpha=0.2)
    ax.fill_between(np.linspace(h0licow-h0licow_err[0] ,h0licow+h0licow_err[1] ) , -10, 10, color='green', alpha=0.2)
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
    

