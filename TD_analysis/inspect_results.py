#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import argparse
import numpy as np
from corner import quantile

import matplotlib.pyplot as plt

# my libs
from TD_analysis.pycs_get_res import *
from Utils.tools import *


def plt_err(config,savefig_dir):

    error = get_error(config)
    
    knts  = [i["knt"] for i in error]
    mltp  = [i["mltype"] for i in error]
    mlcnf = [i["mlconfig"] for i in error]
    
    _       = [i["error"].create_error() for i in error]
    Err_tot = [i["error"].tot for i in error]
    
    # we want to plot: 
    # error vs - knts  - mltype  -mllcs  -mldof
    
    # err vs knts mltype and mllcs
    fig,axes = plt.subplots(3,figsize=(12,12)) 
    lets     = ["AB","AC","BC"]
    comb     = []
    dof_comb = []
    
    for i in range(len(Err_tot)):
        for j in range(len(Err_tot[0])):
            ax = axes[j]
            if mltp[i]=="polyml":
                col  = "r"
                desc = "PolyMl"
                cmb  = "P"
            else:
                col  = "lime"
                desc = "SplMl"
                cmb  = "S"        
            marker="1" #no ml
            dof_mrk = "P"
            alpha=.1
            if mlcnf[i][0]=="AB":
                marker="o"
                cmb+="ab"
            elif mlcnf[i][0]=="AC":
                marker="v"
                cmb+="ac"
            elif mlcnf[i][0]=="BC":
                marker="^"
                cmb+="bc"
            elif mlcnf[i][0]=="ABC":
                marker="s"
                cmb+="abc"
                alpha=1
            if cmb not in comb:
                desc+=" "+str(mlcnf[i][0])
                comb.append(cmb)
            else:
                desc=None 
            ax.set_title(lets[j]+ " Error vs knots")
            ax.set_xlabel("N* knots []")
            ax.set_ylabel("Tot uncertainty [d]")
            ax.scatter(knts[i],Err_tot[i][j],c=col,marker=marker,label=desc,alpha=alpha)
            if mlcnf[i][1]==2:
                dof_mrk = "+"
                dof_desc = "ML dof = 2"
                dof_cmb = "2"
            elif mlcnf[i][1]==3:
                dof_mrk = "x"
                dof_desc = "ML dof = 3"
                dof_cmb = "3"
            else:
                dof_desc ="ML dof = 0"
                dof_cmb  = "0"
            if dof_cmb not in dof_comb:
                dof_comb.append(dof_cmb)
            else:
                dof_desc=None 
            ax.scatter(knts[i],Err_tot[i][j],c="k",marker=dof_mrk,label=dof_desc)
                
            if desc is not None or dof_desc is not None:
                ax.legend()
    plt.tight_layout()
    plt.savefig(savefig_dir+"/err_vs_knt.pdf")
    print("Plotted "+savefig_dir+"/err_vs_knt.pdf")
    plt.close()

def plt_err_tot(config,savefig_dir):
    
    error = get_error(config)
    
    knts  = [i["knt"] for i in error]
    mltp  = [i["mltype"] for i in error]
    mlcnf = [i["mlconfig"] for i in error]
    
    _       = [i["error"].create_error() for i in error]
    Err_tot = [i["error"].tot for i in error]
    lets     = ["AB","AC","BC"]
    
    fig,axes = plt.subplots(1,figsize=(12,6)) 
    comb     = []
    dof_comb = []
    mean_err = []
    mean_err_part = []
    for i in range(len(Err_tot)):
        ax = axes
        if mltp[i]=="polyml":
            col  = "r"
            desc = "PolyMl"
            cmb  = "P"
        else:
            col  = "lime"
            desc = "SplMl"
            cmb  = "S"
        marker="1" #no ml
        dof_mrk = "P"
        alpha = .1
        if mlcnf[i][0]=="AB":
            marker="o"
            cmb+="ab"
        elif mlcnf[i][0]=="AC":
            marker="v"
            cmb+="ac"
        elif mlcnf[i][0]=="BC":
            marker="^"
            cmb+="bc"
        elif mlcnf[i][0]=="ABC":
            marker="s"
            cmb+="abc"
            alpha=1
        if cmb not in comb:
            desc+=" "+str(mlcnf[i][0])
            comb.append(cmb)
        else:
            desc=None
        if mlcnf[i][1]==2:
            dof_mrk = "+"
            dof_desc = "ML dof = 2"
            dof_cmb = "2"
        elif mlcnf[i][1]==3:
            dof_mrk = "x"
            dof_desc = "ML dof = 3"
            dof_cmb = "3"
        else:
            dof_desc ="ML dof = 0"
            dof_cmb  = "0"
        if dof_cmb not in dof_comb:
            dof_comb.append(dof_cmb)
        else:
            dof_desc=None 
        Err_tot_av = np.mean(Err_tot[i])
        #if alpha>.9:
        ax.scatter(knts[i],Err_tot_av,c=col,marker=marker,label=desc) 
        ax.scatter(knts[i],Err_tot_av,c="k",marker=dof_mrk,label=dof_desc,alpha=1)
        
        mean_err.append([knts[i],Err_tot_av])
        if alpha>.9:
            mean_err_part.append([knts[i],Err_tot_av])
    knts_once = []
    for k in np.sort(knts):
        if k not in knts_once:
            knts_once.append(k)
    mean_err_ = [[k,[]] for k in knts_once]
    mean_err_part_ = [[k,[]] for k in knts_once]

    for mr in mean_err:
        knt,merr = mr
        for i_kn in range(len(mean_err_)):
            if knt==mean_err_[i_kn][0]:
                mean_err_[i_kn][1].append(merr)
    mean_err_ = np.transpose(mean_err_)[1].tolist()
    mean_err = [np.mean(mean_err_i) for mean_err_i in mean_err_]
    for mr in mean_err_part:
        knt,merr = mr
        for i_kn in range(len(mean_err_part_)):
            if knt==mean_err_part_[i_kn][0]:
                mean_err_part_[i_kn][1].append(merr)
    mean_err_part_ = np.transpose(mean_err_part_)[1].tolist()
    mean_err_part =  [np.mean(mean_err_part_i).tolist() for mean_err_part_i in mean_err_part_]# np.mean(mean_err_part_,axis=0)
    ax.plot(knts_once,mean_err,alpha=.5,c="gold",label="average")
    ax.plot(knts_once,mean_err_part,alpha=.5,c="red",label="average of selected")
    ax.set_title("Tot Error vs knots")
    ax.set_xlabel("N* knots []")
    ax.set_ylabel("Tot uncertainty [d]")
    ax.legend()
    plt.tight_layout()
    plt.savefig(savefig_dir+"/tot_err_vs_knt.pdf")
    plt.close()


    """
    ###
    # err vs mldof and mltype
    
    fig,axes = plt.subplots(3,figsize=(12,12)) 
    comb     = []
    for i in range(len(Err_tot)):
        for j in range(len(Err_tot[0])):
            ax = axes[j]
            if mltp[i]=="polyml":
                col = "r"
                desc = "PolyMl"
                cmb  = "P"
            else:
                col = "b"
                desc = "SplMl"
                cmb  = "S"
                
            if cmb in comb:
                desc = None
            else:
                comb.append(cmb)
            ax.set_title(lets[j]+ " Error vs DoF")
            ax.set_ylabel("DoF of the ML correction")    
            ax.scatter(Err_tot[i][j],np.transpose(mlcnf)[1][i],c=col,label=desc)
            ax.set_xlabel("Tot uncertainty [d]")
            if desc is not None:
                ax.legend()
    plt.tight_layout()

    plt.savefig(savefig_dir+"/err_vs_dof.png")
    plt.close()
    
    print("Plots produced in "+savefig_dir)




"""
def plt_intr_err(config,savefig_dir):
        fig,axes = plt.subplots(3,figsize=(12,12)) 
        fig.suptitle("Intrinsic scatter from the $\Delta t$ analysis")
        lets     = ["AB","AC","BC"]
        savefig_path = config.analysis_directory
        legends  = []
        for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
            saveknt_path = savefig_path+str("/analysis_kn"+str(knt))
            mkdir(savefig_path)  
            for mlt_i, mltype_i in enumerate(config.mltype): #ml type
                for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                    saveml_path = saveknt_path+"/"+config.get_savemlpath(mltype_i,ml_config)#saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)            
                    #if not check_success_analysis(get_analdir(saveml_path)):
                    #    print("Analysis from "+get_analdir(saveml_path)+" was not sufficiently precise. Ignored .") 
                    #    continue    
                                    
                    with open(str(saveml_path)+"/td.data", 'rb') as f:
                        timedelays=pickle.load(f)
                    distr = np.transpose(timedelays)
                    
                    if mltype_i=="polyml":
                        col  = "r"
                        desc = "PolyMl"
                    else:
                        col  = "lime"
                        desc = "SplMl"
                    marker="1" #no ml
                    dof_mrk = "P"
                    if ml_config[0]=="AB":
                        marker="o"
                        alpha=.1    
                    elif ml_config[0]=="AC":
                        marker="v"
                        alpha=.1
                    elif ml_config[0]=="BC":
                        marker="^"
                        alpha=.1
                    elif ml_config[0]=="ABC":
                        marker="s"
                        alpha=1
                    else:
                        desc=None
                    if ml_config[1]==2:
                        dof_mrk = "+"
                        dof_desc = "ML dof = 2"
                    if ml_config[1]==3:
                        dof_mrk = "x"
                        dof_desc = "ML dof = 3"
                    else:
                        dof_mrk = "P"
                        dof_desc=None 
                    
                    for i in range(len(distr)): #ABCD          
                        ddt = distr[i]
                        #td_i = np.std(ddt)
                        std_i = np.std(ddt)
                        vl,vc,vr = quantile(ddt,q=[0.16,0.5,0.84])
                        err_min = vc-vl
                        err_max = vr-vc
                        mean_err =(err_min+err_max)/2.

                        #axes[i].scatter(knt,std_i,c=col,marker=marker,label=desc) 
                        axes[i].scatter(knt,mean_err,c=col,marker="x",label=desc+" mean error") 
                        axes[i].scatter(knt,std_i,c="k",marker=dof_mrk,label=dof_desc)
                        if a==0 and mlt_i==0:
                            axes[i].set_xlabel("Dt knts [d]")
        for i,ax in enumerate(axes):
            ax.set_title(lets[i])
            ax.set_ylabel("Scatter from time delay [d]")
        plt.tight_layout()
        plt.savefig(savefig_dir+"/scatter_dt_res.pdf")
        print("Plotted: "+savefig_dir+"/scatter_dt_res.pdf")


#"J1433_forcen"
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser        = argparse.ArgumentParser(description="plot error wrt config",
                                   formatter_class=argparse.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process (Euler, SMARTS, ... )"

    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-err','--error_analysis',help="Plot the error per knot, per lcs",
                        dest="err", 
                        default=False,action="store_true")
    parser.add_argument('-err_tot','--error_tot_analysis',help="Plot the total error per knot, per lcs",
                        dest="err_tot", 
                        default=False,action="store_true")
    parser.add_argument('-intr_err','--intrinsic_error_analysis',help="Plot the scatter of the dt analysis (intrinsic error) per knot, per lcs",
                        dest="intr_err", 
                        default=False,action="store_true")
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname
    err      = args.err
    err_tot  = args.err_tot
    intr_err = args.intr_err   
    
    config = lensname+"_"+dataname
    config = get_config(config)
    savefig_dir = config.combined_directory +'/figure/marginalisation_plots/'
    if err:
        plt_err(config,savefig_dir)
    
    if err_tot:
        plt_err_tot(config,savefig_dir)
        
    if intr_err:
        plt_intr_err(config,savefig_dir)
        """
        savefig_path = config.analysis_directory
        ind_col = 0 
        colors = config.colors
        legend_elements = []
        knts = []
        stdv = []
        for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
            saveknt_path = savefig_path+str("/analysis_kn"+str(knt))
            mkdir(savefig_path)  
            stdv_mean = []
            knts.append(knt)
            for mlt_i, mltype_i in enumerate(config.mltype): #ml type
                for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                    mllist_name,mlfp = ml_config 
                    if mltype_i=="polyml":
                        mlfp_str="_mlfp_"+str(mlfp)
                    else:
                        if config.forcen:
                            mlfp_str="_nmlspl_"+str(mlfp)
                        else:
                            mlfp_str="_knstml_"+str(mlfp)
                            
                    saveml_path = saveknt_path+str("/ml"+mltype_i[:-2]+mllist_name+mlfp_str)

                    with open(str(saveml_path)+"/td.data", 'rb') as f:
                        timedelays=pickle.load(f)
                    distr = np.transpose(timedelays)
                    stdv_ABC = []
                    for i in range(len(distr)): #ABCD          
                        ddt = distr[i]
                        std_i = np.std(ddt)
                        stdv_ABC.append(std_i)
                    stdv_mean.append(stdv_ABC)
            stdv.append(np.mean(stdv_mean,axis=0))
        
        fg,ax = plt.subplots(3,figsize=(20,20))
        for i in range(len(ax)):
            ax[i].scatter(knts,np.transpose(stdv)[i])
        plt.savefig(savefig_dir+"/scatter_dt_res.pdf")
        """

    success(sys.argv[0])
