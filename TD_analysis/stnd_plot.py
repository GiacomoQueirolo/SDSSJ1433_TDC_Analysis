#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# In[ ]:


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
        
def delayplot(groups,savename,colors,refgroup=None,selected_groups_indexes=[]):
    print("Note: we are considering only lc A,B and C")
    fg,ax= plt.subplots(2,2,figsize=(15,15))  

    prm_df = [r"$\Delta t$ "+lbl+" [d]" for lbl in groups[0].labels]

    def fmt_new(y,pos):
        return ""

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

        fntsz = 10
        for j in range(len(groups)):
            G = groups[j]
            
            #td_err_i = [[G.errors_down[i]],[G.errors_up[i]]] method A
            td_err_i = G.tot_error[i]
            td_res_i = G.results[i]
            if j >= len(colors):
                colj=colors[int(j%len(colors))]
            else:
                colj=colors[j]
            ax_i.errorbar(td_res_i,y_rnd[j],yerr=None,xerr=td_err_i,fmt=colj,capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[j], c=colj)
            sel_txt= ""
            if j in selected_groups_indexes:
                cir = plt.Circle( (td_res_i,y_rnd[j]), 18/plt.rcParams['figure.dpi'], color=colj,fill=False)
                ax_i.set_aspect('equal', adjustable='datalim')
                ax_i.add_patch(cir)
                sel_txt=" (S)"
            str_res = print_res_w_err(td_res_i,td_err_i)+sel_txt
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[j]+.2,s=str_res,c=colj,fontsize=fntsz)
        ax_i.set_xlabel("$\Delta t$ [d]")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

        if refgroup is not None:
            #td_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            td_err_i = refgroup.tot_error[i]
            td_res_i = refgroup.results[i] 
            ax_i.errorbar(td_res_i,y_rnd[-1],yerr=None,xerr=td_err_i,fmt="k",capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[-1], c="k")
            str_res = print_res_w_err(td_res_i,td_err_i)
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[-1]+.2,s=str_res,c="k",fontsize=fntsz)
            ax_i.fill_between([td_res_i-td_err_i,td_err_i+td_res_i], -100, 100, color='grey', alpha=0.2) 
    ax_del = ax[0][1]
    ax_del.axis("off")
    ln_lgnd = []
    for j in range(len(groups)):
        G = groups[j]        
        if j >= len(colors):
            colj=colors[int(j%len(colors))]
        else:
            colj=colors[j]
        y_lbl = G.name.replace("_"," ")
        ln, = ax_del.plot(1,1,c=colj,label=y_lbl)
        ln_lgnd.append(ln)
    if refgroup is not None:
        if len(str(refgroup.name))>30:
            ln, = ax_del.plot(1,1,c="k",label="Multiple groups combined")
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
    
    plt.savefig(savename)
    print("Saving plot "+savename)
    


# In[ ]:


def dmagplot(groups,savename,colors,refgroup=None):
    print("Note: we are considering only lc A,B and C")
    fg,ax= plt.subplots(2,2,figsize=(13,13))  

    prm_df = [r"$\Delta mag$ "+lbl+" [d]" for lbl in groups[0].labels]

    def fmt_new(y,pos):
        return ""

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

        fntsz = 10
        for j in range(len(groups)):
            G = groups[j]
            
            #dmag_err_i = [[G.errors_down[i]],[G.errors_up[i]]] method A
            dmag_err_i = G.tot_error[i]
            dmag_res_i = G.results[i]
            if j >= len(colors):
                colj=colors[int(j%len(colors))]
            else:
                colj=colors[j]
            ax_i.errorbar(dmag_res_i,y_rnd[j],yerr=None,xerr=dmag_err_i,fmt=colj,capsize=4)
            ax_i.scatter(dmag_res_i ,y_rnd[j], c=colj)
            str_res = print_res_w_err(dmag_res_i,dmag_err_i)
            ax_i.text(dmag_res_i,y_rnd[j]+.2,s=str_res,c=colj,fontsize=fntsz)
        
        ax_i.set_xlabel("$\Delta mag$ []")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

        if refgroup is not None:
            #dmag_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            dmag_err_i = refgroup.tot_error[i]
            dmag_res_i = refgroup.results[i] 
            ax_i.errorbar(dmag_res_i,y_rnd[-1],yerr=None,xerr=dmag_err_i,fmt="k",capsize=4)
            ax_i.scatter(dmag_res_i ,y_rnd[-1], c="k")
            str_res = print_res_w_err(dmag_res_i,dmag_err_i)
            ax_i.text(dmag_res_i,y_rnd[-1]+.2,s=str_res,c="k",fontsize=fntsz)
            ax_i.fill_between([dmag_res_i-dmag_err_i,dmag_err_i+dmag_res_i], -100, 100, color='grey', alpha=0.2) 

    ax_del = ax[0][1]
    ax_del.axis("off")
    ln_lgnd = []
    for j in range(len(groups)):
        G = groups[j]        
        if j >= len(colors):
            colj=colors[int(j%len(colors))]
        else:
            colj=colors[j]
        y_lbl = G.name.replace("_"," ")
        ln, = ax_del.plot(1,1,c=colj,label=y_lbl)
        ln_lgnd.append(ln)
    if refgroup is not None:
        if len(str(refgroup.name))>30:
            ln, = ax_del.plot(1,1,c="k",label="Multiple groups combined")
        else:
            ln, = ax_del.plot(1,1,c="k",label=refgroup.name)
        ln_lgnd.append(ln)
        
    if len(ln_lgnd)<20:
        ax_del.legend()   
    else:
        set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
        ax_del.add_artist(set_lgnd)
        ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")
    
    plt.savefig(savename)
    print("Saving plot "+savename)

