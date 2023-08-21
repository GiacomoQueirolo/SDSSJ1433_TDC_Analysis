#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np

from H0.tools_H0 import *
from Utils.statistical_tools import quantiles


# In[2]:


font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}
matplotlib.rc('font', **font)
matplotlib.rc('figure',**{'figsize':(12,9)})
plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)


# In[4]:


def plot_H0(H0,PH0,figname=False,title=r"$P(H_0)$",return_plot=False,add_mode=True,confront=True,figsize=(12,12),color="b",ax=None):
    """
    # Plot H0 posterior compared to the two main results and 1-sigma regions of the litterature (H0LiCOW and Planck)
    # Input:
    # H0: array
    #     sampling of H0
    # PH0: array
    #     posterior of H0
    # figname: string (optional)
    #     name of the saved plot (if to be saved)
    """
    x = H0
    y = PH0
    yh0 = max(PH0)/2
    h0_res,err_min,err_max= quantiles(PH0,H0,return_quantiles=False)
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=figsize) 
    ax.plot(x,y,'+:',color=color)
    ax.errorbar(h0_res,yh0,yerr=None,xerr=[[err_min],[err_max]],fmt=color,capsize=4)
    ax.scatter(h0_res ,yh0, c=color,marker="*",s=100)#,label="$H_0$"
    string_s=str(np.round(h0_res,2))+"$_{-"+str(np.round(err_min,2))+"}^{+"+str(np.round(err_max,2))+"}$"
    ax.text(h0_res-.5*len(str(np.round(h0_res,3))),yh0*1.05,s=string_s,c="k",fontsize=16)
    ax.set_title(title)
    ax.set_ylabel("Probability")
    ax.set_xlabel(r"$H_0\ [km/s/Mpc]$")
    ax.set_ylim(0,max(y)*1.1)
    if add_mode:
        ind_mode = np.where(PH0==np.max(PH0))[0][0]
        mode     = H0[ind_mode]
        yh0_mode = max(PH0)/3
        ax.scatter (mode,yh0_mode,c="r",marker="*",s=100)
        ax.errorbar(mode,yh0_mode,yerr=None,xerr=[[err_min],[err_max]],fmt="r",capsize=4)
        ax.errorbar(None,None,fmt="k",label="Median")
        ax.errorbar(None,None,fmt="r",label="Mode")
        string_mode = str(np.round(mode,2))+"$_{-"+str(np.round(err_min,2))+"}^{+"+str(np.round(err_max,2))+"}$"
        ax.text(mode-.5*len(str((np.round(mode,2)))),yh0_mode*1.05,string_mode)    
        ax.legend()
    if confront:
        ax.axvline(h0planck.H0,label="Planck", c="r")
        ax.axvline(h0licow.H0, label="H0LiCOW",c="g")
        ax.fill_between(np.linspace(h0planck.H0-h0planck.sigma_min, h0planck.H0+ h0planck.sigma_max), -10, 10, color='r', alpha=0.2)
        ax.fill_between(np.linspace(h0licow.H0-h0licow.sigma_min,   h0licow.H0 + h0licow.sigma_max),  -10, 10, color='g', alpha=0.2)
        ax.legend()  
    if return_plot:
        return ax
    if type(figname) is str:
        plt.savefig(figname, transparent=True)
    else:
        plt.show()

