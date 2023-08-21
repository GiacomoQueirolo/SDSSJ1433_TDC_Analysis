#!/usr/bin/env python
# coding: utf-8
#I want to analyse the stacked MCMC run to re-write the results parameter consistently with the MCMC results

import sys
import pickle
import argparse
from corner import quantile

from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.get_res_multifilter import get_mcmc_smpl_mltf,get_mcmc_prm_mltf

 
@check_mltf_setting
def rewrite_read_results_mltf(multifilter_setting,cut_mcmc=0):
    #MCMC sample
    samples_mcmc = get_mcmc_smpl_mltf(multifilter_setting)[cut_mcmc:]
    #parameters' name
    param_mcmc   = get_mcmc_prm_mltf(multifilter_setting)[cut_mcmc:]
    
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

    with open(f"{multifilter_setting.savemcmc_path}/read_results_updated.data", 'wb') as file:
        pickle.dump(kwargs_results_updated, file)
 

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Produces the read_results_updated.data with the correct values relative to the MCMC run - MULTIFILTER_VERSION")
    parser.add_argument('MULTIFILTER_SETTING',nargs="1",default=[],help="Multifilter_setting file to consider")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    args = parser.parse_args()
    multifilter_setting = get_multifilter_setting_module(args.MULTIFILTER_SETTING)
    cut_mcmc            = args.cut_mcmc

    present_program(sys.argv[0])
    rewrite_read_results_mltf(multifilter_setting,cut_mcmc)
    success(sys.argv[0])

