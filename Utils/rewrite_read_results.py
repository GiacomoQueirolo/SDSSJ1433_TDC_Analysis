#!/usr/bin/env python
# coding: utf-8
#I want to analyse the stacked MCMC run to re-write the results parameter consistently with the MCMC results
# In[ ]:


import sys
import pickle
import argparse
import numpy as np
from corner import quantile

from Utils.tools import *
from Utils.get_res import load_whatever


# In[ ]:


def rewrite_read_results(setting,cut_mcmc=0,backup_path="backup_results"):
    #MCMC sample
    samples_mcmc = get_mcmc_smpl(setting,backup_path)[cut_mcmc:]
    #parameters' name
    param_mcmc   = get_mcmc_prm(setting,backup_path)
    
    kwargs_results_updated={} 
    n_ra  = 0
    n_dec = 0  

    for i in range(len(param_mcmc)):
        val = float(quantile(samples_mcmc[:,i],q=.5))
        if param_mcmc[i]!="ra_image" and param_mcmc[i]!="dec_image":
            kwargs_results_updated[param_mcmc[i]]=val 
        elif param_mcmc[i]=="ra_image":
            kwargs_results_updated["ra_image_"+str(n_ra)]=val          
            n_ra+=1
        else:
            kwargs_results_updated["dec_image_"+str(n_dec)]=val          
            n_dec+=1

    with open(get_savemcmcpath(setting_name,backup_path)+'/read_results_updated.data', 'wb') as file:
        pickle.dump(kwargs_results_updated, file)


# In[ ]:


if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Produces the read_results_updated.data with the correct values relative to the MCMC run")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    args = parser.parse_args()
    settings = args.SETTING_FILE
    cut_mcmc = args.cut_mcmc
    present_program(sys.argv[0])
    for setting in settings:
        rewrite_read_results(setting,cut_mcmc,backup_path)
    success(sys.argv[0])

