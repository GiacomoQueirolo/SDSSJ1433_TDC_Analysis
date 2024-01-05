#!/usr/bin/env python
# coding: utf-8

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

comb_label = "Combined result"

matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['text.usetex'] = True 

def print_res_w_err(res,err,outstr=True):
    exp=int(str("{:.1e}".format(err)).split("e")[-1])
    if exp>-2 and exp<1:
        exp=-1
    str_res = str(np.round(res,abs(exp)))+"$\pm$"+str(np.round(err,abs(exp)))
    # to implement in case of very large/small numbers to use the scientific notation -> not the case here
    if outstr:
        return str_res
    else:
        print(str_res)
        
def rename_label(group_name):
    readable_label = r""
    group_name = str(group_name).replace("_"," ")
    kntstp = group_name.split("ks")[1].split(" ")[0]
    readable_label = "Knots. = "+kntstp
    if "ml" in group_name:
        if "polyml" in group_name:
            deg = group_name.split("deg ")[1][0]
            readable_label+=", "+deg+r"$^o$ Polynomial ML"
        elif "splml" in group_name:
            nmspl = group_name.split("nmlspl ")[1][0]
            readable_label+=", "+nmspl+"-piece Spline ML"
        letts=["A","B","C","D"]
        if all([True if let in group_name else False for let in letts]):
            readable_label+=r" $\forall$ images"
        else:
            readable_label+=" applied to "
            for let in letts:
                if let in group_name:
                    readable_label+=let
                    if let!="D":
                        readable_label+=", "
    else:
        readable_label+=" without ML correction"
    return readable_label


def fmt_new(y,pos):
        return ""

def delayplot(groups,savename=None,refgroup=None,selected_groups_indexes=[],readable_labels=True):
    fntsz = 30
    plt.rcParams['xtick.labelsize'] = fntsz
    plt.rcParams['ytick.labelsize'] = fntsz 
    plt.rcParams['font.size'] = fntsz
    plt.rc('axes', labelsize=fntsz)     # fontsize of the x and y labels
    plt.rc('font', size=fntsz)          # controls default text sizes
    plt.rc('axes', titlesize=fntsz)     # fontsize of the axes title
    plt.rc('axes', labelsize=fntsz)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fntsz)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fntsz)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fntsz)    # legend fontsize
    fg,ax = plt.subplots(2,2,figsize=(16,16))  

    prm_df = [r"$\Delta t$ "+lbl+" [d]" for lbl in groups[0].labels]

    
    for i in range(len(prm_df)):

        if i==0:
            ax_i=ax[0][0]
        if i==1:
            ax_i=ax[1][0]
        if i==2:
            ax_i=ax[1][1]
        ax_i.yaxis.set_major_formatter(FuncFormatter(fmt_new))
        
        y_rnd = [0.5 + j for j in range(len(groups),0,-1)]
        if refgroup is not None:
            y_rnd.append(y_rnd[-1]-1)

        
        for j in range(len(groups)):
            G = groups[j]
            
            #td_err_i = [[G.errors_down[i]],[G.errors_up[i]]] method A
            td_err_i = G.tot_error[i]
            td_res_i = G.results[i]
            colj     = G.color
            ax_i.errorbar(td_res_i,y_rnd[j],yerr=None,xerr=td_err_i,fmt=colj,capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[j], c=colj)
            sel_txt= ""
            if refgroup:
                if G.name in refgroup.combined_names:
                    cir = plt.Circle( (td_res_i,y_rnd[j]), 18/plt.rcParams['figure.dpi'], color=colj,fill=False)
                    if i!=2:
                        ax_i.set_aspect('equal', adjustable='datalim')
                    ax_i.add_patch(cir)
                    sel_txt=" (S)"
            str_res = print_res_w_err(td_res_i,td_err_i)+sel_txt
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[j]+.2,s=str_res,c=colj,fontsize=fntsz)

        if refgroup is not None:
            #td_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            td_err_i = refgroup.tot_error[i]
            td_res_i = refgroup.results[i] 
            ax_i.errorbar(td_res_i,y_rnd[-1],yerr=None,xerr=td_err_i,fmt="k",capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[-1], c="k")
            str_res = print_res_w_err(td_res_i,td_err_i)
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[-1]+.2,s=str_res,c="k",fontsize=fntsz)
            ax_i.fill_between([td_res_i-td_err_i,td_err_i+td_res_i], -100, 100, color='grey', alpha=0.2) 
        ax_i.set_xlabel("$\Delta t$ [d]")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)
    ax_del = ax[0][1]
    ax_del.axis("off")
    ln_lgnd = []
    for j in range(len(groups)):
        G = groups[j]        
        colj=G.color
        y_lbl = G.name.replace("_"," ")
        if readable_labels:
            y_lbl = rename_label(G.name)
        ln, = ax_del.plot(1,1,c=colj,label=y_lbl)
        ln_lgnd.append(ln)
    if refgroup is not None:
        if len(str(refgroup.name))>30:
            ln, = ax_del.plot(1,1,c="k",label=comb_label)
        else:
            ln, = ax_del.plot(1,1,c="k",label=refgroup.name)
        ln_lgnd.append(ln)
    if len(selected_groups_indexes)!=0:
        ln, = ax_del.plot(1,1,c="w",label="(S): Selected groups for combined result")
        ln_lgnd.append(ln)
        
    if len(ln_lgnd)<20:
        ax_del.legend()   
    else:
        set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
        ax_del.add_artist(set_lgnd)
        ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")
    if savename:
        plt.tight_layout()
        plt.savefig(savename)
        print("Saving plot "+savename)
    else:
        return fg
        
def dmagplot(groups,savename=None,refgroup=None,readable_label=True):
    fntsz = 30
    plt.rcParams['xtick.labelsize'] = fntsz
    plt.rcParams['ytick.labelsize'] = fntsz 
    plt.rcParams['font.size'] = fntsz
    plt.rc('axes', labelsize=fntsz)     # fontsize of the x and y labels
    plt.rc('font', size=fntsz)          # controls default text sizes
    plt.rc('axes', titlesize=fntsz)     # fontsize of the axes title
    plt.rc('axes', labelsize=fntsz)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fntsz)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fntsz)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fntsz)    # legend fontsize
    fg,ax = plt.subplots(2,2,figsize=(17,17))  
    
    prm_df = [r"$\Delta $mag "+lbl+" [ ]" for lbl in groups[0].labels]

    for i in range(len(prm_df)):

        if i==0:
            ax_i=ax[0][0]
        if i==1:
            ax_i=ax[1][0]
        if i==2:
            ax_i=ax[1][1]
        ax_i.yaxis.set_major_formatter(FuncFormatter(fmt_new))
        
        y_rnd = [0.5 + j for j in range(len(groups),0,-1)]
        if refgroup is not None:
            y_rnd.append(y_rnd[-1]-1)

        for j in range(len(groups)):
            G = groups[j]
            
            #dmag_err_i = [[G.errors_down[i]],[G.errors_up[i]]] method A
            dmag_err_i = G.tot_error[i]
            dmag_res_i = G.results[i]
            colj       = G.color
            ax_i.errorbar(dmag_res_i,y_rnd[j],yerr=None,xerr=dmag_err_i,fmt=colj,capsize=4)
            ax_i.scatter(dmag_res_i ,y_rnd[j], c=colj)
            str_res = print_res_w_err(dmag_res_i,dmag_err_i)
            ax_i.text(dmag_res_i,y_rnd[j]+.2,s=str_res,c=colj,fontsize=fntsz)

        if refgroup is not None:
            #dmag_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            dmag_err_i = refgroup.tot_error[i]
            dmag_res_i = refgroup.results[i] 
            ax_i.errorbar(dmag_res_i,y_rnd[-1],yerr=None,xerr=dmag_err_i,fmt="k",capsize=4)
            ax_i.scatter(dmag_res_i ,y_rnd[-1], c="k")
            str_res = print_res_w_err(dmag_res_i,dmag_err_i)
            ax_i.text(dmag_res_i,y_rnd[-1]+.2,s=str_res,c="k",fontsize=fntsz)
            ax_i.fill_between([dmag_res_i-dmag_err_i,dmag_err_i+dmag_res_i], -100, 100, color='grey', alpha=0.2) 
        ax_i.set_xlabel("$\Delta $mag [ ]")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

    ax_del = ax[0][1]
    ax_del.axis("off")
    ln_lgnd = []
    for j in range(len(groups)):
        G     = groups[j]        
        colj  = G.color
        y_lbl = G.name.replace("_"," ")
        if readable_label:
            y_lbl = rename_label(G.name)
        ln, = ax_del.plot(1,1,c=colj,label=y_lbl)
        ln_lgnd.append(ln)
    if refgroup is not None:
        if len(str(refgroup.name))>30:
            ln, = ax_del.plot(1,1,c="k",label=comb_label)
        else:
            ln, = ax_del.plot(1,1,c="k",label=refgroup.name)
        ln_lgnd.append(ln)
        
    if len(ln_lgnd)<20:
        ax_del.legend()   
    else:
        set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
        ax_del.add_artist(set_lgnd)
        ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")
    if savename:
        plt.tight_layout()
        plt.savefig(savename)
        print("Saving plot "+savename)
    else:
        return fg
