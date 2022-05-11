#!/usr/bin/env python
# coding: utf-8

import corner
import pickle
import json
import numpy as np
import argparse
import matplotlib.lines as mlines
from my_lenstronomy.my_kwargs2args import my_kwarg2args

from tools import *
from get_res import * 

from fermat_pot_analysis import gen_mcmc_Df,labels_Df


def MedianVsLogL(setting_name,cut_mcmc):
    setting       = get_setting_module(setting_name).setting()
    savefig_path  = get_savefigpath(setting_name,backup_path)
    savemcmc_path = get_savemcmc_path(setting_name,backup_path)

    ws = check_if_WS(setting_name)
    #MCMC sample
    samples_mcmc = get_mcmc_smpl(setting_name,backup_path)
    #parameters' name
    param_mcmc = get_mcmc_prm(setting_name,backup_path)
    
    labels = []
    for prm in param_mcmc:
        if "image" not in prm:
            nm,udm = setting.str_param(prm)
        else:
            nm = prm
            udm = "[\"]"
        labels.append(nm+" "+udm)
        
    if cut_mcmc>0:
        samples_mcmc = samples_mcmc[cut_mcmc:]        

    # results from median of each param
    res_med =np.mean(samples_mcmc,axis=0)

    # results from best logL
    kwargs_res = get_kwres(setting_name,backup_path)["kwargs_results"]

    res_logL =  my_kwarg2args(kwargs_res,setting_name,cut_mcmc=cut_mcmc)


    #Plot for All Params
    figure = corner.corner(samples_mcmc,labels=labels,show_titles=True)

    corner.overplot_lines(figure, res_med, color="red")
    corner.overplot_points(figure, res_med[None], marker="s", color="red")
    corner.overplot_lines(figure, res_logL, color="blue")
    corner.overplot_points(figure, res_logL[None], marker="s", color="blue")

    line_med = mlines.Line2D([], [], color='red',label='Median')
    line_logL = mlines.Line2D([], [], color='blue',label='Best logL')
    figure.legend(handles=[line_med,line_logL])

    figure.savefig(savefig_path+"/MedianVsLogL.png")
    
    ###################
    # Fermat pot
    labels_Df = ["$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} AD$"]
    
    mcmc_f   = get_mcmc_fermat(setting_name,backup_path) 
    mcmc_DfT = np.transpose(mcmc_f)[1:]-np.transpose(mcmc_f)[0] 
    mcmc_Df  = mcmc_DfT.T
    #D_AB, D_AC, D_AD, meaning Df_i - Df_A
    if cut_mcmc>0:
        mcmc_Df = mcmc_Df[cut_mcmc:]
        
    # results from median of each param
    res_med_Df =np.mean(mcmc_Df,axis=0)

    # results from best logL
    # same as before
    # convert it to fermat pot as a single step mcmc
    

    res_logL_Df = gen_mcmc_Df([res_logL,res_logL],param_mcmc,setting)[0]


    #Plot for All Params
    figure = corner.corner(mcmc_Df,labels=labels_Df,show_titles=True)

    corner.overplot_lines(figure, res_med_Df, color="red")
    corner.overplot_points(figure, res_med_Df[None], marker="s", color="red")
    corner.overplot_lines(figure, res_logL_Df, color="blue")
    corner.overplot_points(figure, np.array(res_logL_Df)[None], marker="s", color="blue")

    line_med  = mlines.Line2D([], [], color='red',label='Median')
    line_logL = mlines.Line2D([], [], color='blue',label='Best logL')
    figure.legend(handles=[line_med,line_logL])

    figure.savefig(savefig_path+"/MedianVsLogL_Df.png")

    
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Produces the corner plot with comparison between read_updated_results.data (Median) and read_results.data (Best LogL)")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings = args.SETTING_FILES
    cut_mcmc =args.cut_mcmc

    backup_path="backup_results"
        
    for sets in settings:
        if len(settings)>1:
            print("Plotting: ",sets)
        MedianVsLogL(sets,cut_mcmc)
    success(sys.argv[0])
