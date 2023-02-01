#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


# In[ ]:


def delayplot(groups,savename,colors,refgroup=None):
    print("Note: we are considering only lc A,B and C")
    fg,ax= plt.subplots(2,2,figsize=(30,30))  

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
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[j]+.2,s=str(np.round(td_res_i,2))+
                 "$\pm$"+str(np.round(td_err_i,2)),c=colj,fontsize=fntsz)
        
        ax_i.set_xlabel("$\Delta t$ [d]")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

        if refgroup is not None:
            #td_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            td_err_i = refgroup.tot_error[i]
            td_res_i = refgroup.results[i] 
            ax_i.errorbar(td_res_i,y_rnd[-1],yerr=None,xerr=td_err_i,fmt="k",capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[-1], c="k")
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[-1]+.2,s=str(np.round(td_res_i,2))+
                      "$\pm$"+str(np.round(td_err_i,2)),c="k",fontsize=fntsz)
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
        ln, = ax_del.plot(1,1,c="k",label=refgroup.name)
        ln_lgnd.append(ln)
    if len(ln_lgnd)<20:
        ax_del.legend()   
    else:
        set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
        ax_del.add_artist(set_lgnd)
        ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")
    
    plt.savefig(savename)


# In[ ]:


def dmagplot(groups,savename,colors,refgroup=None):
    print("Note: we are considering only lc A,B and C")
    fg,ax= plt.subplots(2,2,figsize=(30,30))  

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
            
            #td_err_i = [[G.errors_down[i]],[G.errors_up[i]]] method A
            td_err_i = G.tot_error[i]
            td_res_i = G.results[i]
            if j >= len(colors):
                colj=colors[int(j%len(colors))]
            else:
                colj=colors[j]
            ax_i.errorbar(td_res_i,y_rnd[j],yerr=None,xerr=td_err_i,fmt=colj,capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[j], c=colj)
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[j]+.2,s=str(np.round(td_res_i,2))+
                 "$\pm$"+str(np.round(td_err_i,2)),c=colj,fontsize=fntsz)
        
        ax_i.set_xlabel("$\Delta mag$ []")
        ax_i.set_title(prm_df[i])
        ax_i.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

        if refgroup is not None:
            #td_err_i = [[refgroup.errors_down[i]],[refgroup.errors_up[i]]] method A
            td_err_i = refgroup.tot_error[i]
            td_res_i = refgroup.results[i] 
            ax_i.errorbar(td_res_i,y_rnd[-1],yerr=None,xerr=td_err_i,fmt="k",capsize=4)
            ax_i.scatter(td_res_i ,y_rnd[-1], c="k")
            ax_i.text(td_res_i-.05*len(str(np.round(td_res_i,3))),y_rnd[-1]+.2,s=str(np.round(td_res_i,2))+
                      "$\pm$"+str(np.round(td_err_i,2)),c="k",fontsize=fntsz)
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
        ln, = ax_del.plot(1,1,c="k",label=refgroup.name)
        ln_lgnd.append(ln)
    if len(ln_lgnd)<20:
        ax_del.legend()   
    else:
        set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
        ax_del.add_artist(set_lgnd)
        ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")
    
    plt.savefig(savename)

