#!/usr/bin/env python
# coding: utf-8

# In[1]:


import json
from Utils.get_res import *
from Utils.tools import *
import argparse
import corner
import matplotlib.pyplot as plt
import numpy as np
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Stack the mag ratio mcmc chains, obtain a single result and plot it")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()

    settings  = args.SETTING_FILES
    save_plot = "backup_results/PDF_superposition_II/mu_stacking/"
    mkdir(save_plot)

    ####################
    present_program(sys.argv[0])
    ####################

    mag_stack = []
    for st in settings:
        mag_stack.append(get_mcmc_mag(st))
    mag_stack = np.vstack(mag_stack).T # shape = 3, n* steps
    labels = ["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_D$/$\mu_A$"]
    for mg_i in range(len(mag_stack)):

        vmin,res,vmax = corner.quantile(mag_stack[mg_i],q=[0.16,0.5,0.84])
        errmin = res-vmin
        errmax = vmax-res
        print("res =",res,"errmin =",errmin,"errmax =",errmax)
        plt.figure(figsize=(12,9))    
        sigma = np.mean([errmin,errmax])        
        count, bins, _ = plt.hist(mag_stack[mg_i], density=True,stacked=True,label=labels[mg_i]+" :"+str(round(res,2))+r"$\pm$"+str(round(sigma,3)) )

        plt.errorbar(res,max(count)/2.,yerr=None,xerr=[[errmin],[errmax]],fmt="+")
        plt.scatter(res,max(count)/2.,marker=".")

        plt.title(r"Stacking of "+labels[mg_i])
        plt.xlabel(labels[mg_i]+" []")
        plt.legend(loc="upper right")
        param = "mu_"+labels[mg_i].replace("\\","").replace("/","").replace("$","").replace("mu_","")[::-1]
        plt.savefig(str(save_plot)+"/"+param+".png")

    success(sys.argv[0])
