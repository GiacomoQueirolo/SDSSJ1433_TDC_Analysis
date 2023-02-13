#!/usr/bin/env python
# coding: utf-8

import numpy as np
import argparse
from Custom_Model.my_lenstronomy.my_kwargs2args import my_kwarg2args

from Utils.tools import *
from Utils.get_res import * 

from Posterior_analysis.fermat_pot_analysis import gen_mcmc_Df,labels_Df


def plot_Dfs(Dfs,names,figsize=(12,12),nticks=6,colors=["r","b","g","y","k"]):
    
    ticks  = []
    for i in range(len(Dfi[0])): # for each dim
        mx = max([Dfi[i] for Dfi in Dfs]) # find max
        mn = min([Dfi[i] for Dfi in Dfs]) # and min
        dm = mx-mn
        lft = mn-0.01*dm
        rgt = mx+0.01*dm
        stp = (rgt - lft)/nticks
        ticks_i = []
        for k in range(nticks):
            ticks_i.append(np.round(lft + k*stp,2))
        ticks.append(ticks_i)
        
    fg,ax= plt.subplots(2,2,figsize=figsize) 
    for k,Df in enumerate(Dfs):
        col= colors[k]
        print("Df ",names[k])
        for i in range(2):
            for j in range(2):
                ax_ij = ax[i][j]
                if [i,j]!= [0,1]:
                    x = j
                    y = i+1
                    """if [i,j] == [0,0]: #AC vs AB
                        x,y = 0,1
                    elif [i,j]==[1,0]: #BC vs AB
                        x,y = 0,2 
                    elif [i,j]==[1,1]:
                        x,y = 1,2"""
                    if i!=j:
                        ax_ij.set_xlabel(labels_Df[x])
                        ax_ij.set_ylabel(labels_Df[y])
                        ax_ij.set_xticks(ticks[x])
                    elif i==0:
                        ax_ij.set_ylabel(labels_Df[y])
                        ax_ij.set_yticks(ticks[y])
                    elif j==1:                    
                        ax_ij.set_xticks(ticks[x])
                        ax_ij.set_xlabel(labels_Df[x])

                    ax_ij.axvline(Df[x],c=col)
                    ax_ij.axhline(Df[y],c=col) #label=names[k]

    
    axdel=ax[0][1]
    for i,set_name in enumerate(names):
        axdel.scatter([1],[1],label=set_name,color=colors[i])
    axdel.scatter([1],[1], color="w")
    axdel.legend()
    axdel.axis("off")  
    fg.suptitle(r"Compare $\Delta \phi$ from different filters", fontsize=70)  
            
    plt.tight_layout()
    return plt

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Compare the results from the PSO for the different filters")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the plot")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings = args.SETTING_FILES
    dir_name = args.dir_name

    backup_path = "backup_results"
    names  = [strip_setting_name(st) for st in settings]
    st_pos = [find_setting_path(st) for st in settings]
    
    
    savefig_path = create_dir_name(settings,save_dir="PDF_superposition_II",dir_name=dir_name,backup_path=backup_path)
    
    Dfs = []            
    for i,sets in enumerate(settings):
        if len(settings)>1:
            print("Considering: ",sets)
        # results from best logL
        kwargs_res  = get_kwres(sets,backup_path)["kwargs_results"]
        param_mcmc  = get_mcmc_prm(sets,backup_path)
        res_logL    = my_kwarg2args(kwargs_res,sets)
        setting     = get_setting_module(sets).setting()
        res_logL_Df = gen_mcmc_Df([res_logL,res_logL],param_mcmc,setting)[0] #D_AB, D_AC, D_AD
        Dfs.append(res_logL_Df)
    plt = plot_Dfs(Dfs,names=names)
    plt.title(r"Compare $\Delta \phi$ from different filters")
    plt.savefig(savefig_path+"/PSO_results_compare.png") 
    print("Plot saved as "+savefig_path+"/PSO_results_compare.png")


    success(sys.argv[0])
