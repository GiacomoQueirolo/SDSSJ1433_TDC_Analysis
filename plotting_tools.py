#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt
from statistical_tools import marginalise_prob
import numpy as np

base_colors = ["r","b","g","y","k","m","c","darkorange","darkviolet","lawngreen","violet"] 

def plot_probability3D(prob3D,bins,labels=None,udm=None,figsize=(12,12),contour_levels=30,colors=base_colors,alpha=1):
    fg,ax= plt.subplots(3,3,figsize=figsize) 
    bin_coord = [[(bins[j][i]+bins[j][i+1])/2. for i in range(0,len(bins[j])-1)]for j in range(len(bins))]
    bin_densities = marginalise_prob(prob3D,bins) 
    ticks = []
    for i in range(len(bins)):
        rng = bins[i][-1]- bins[i][0]
        stp = bins[i][1] - bins[i][0]
        ticks_i=[]
        for j in range(len(bins[i])):
            ticks_i.append(np.round(bins[i][0]+j*stp,1))
        ticks.append(ticks_i)
        
    for i in range(3):
        for j in range(3):
            ax_ij = ax[i][j]
            if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                if i==j:
                    if labels is not None:
                        if i==0:
                            ax_ij.set_ylabel(labels[i])
                        elif i==2:
                            ax_ij.set_xlabel(labels[i])
                    cnt = bin_densities[i]
                    bn  = bins[i]
                    for k in range(len(cnt)):
                        ax_ij.fill_between([bn[k],bn[k+1]],[0],cnt[k],color=colors[i],alpha=alpha)
                    ax_ij.set_xticks(ticks[j])
                    ax_ij.set_ylim(0,max(cnt)+.2*max(cnt))
                    mx = np.where(cnt==max(cnt))[0][0]
                    #res = str(np.round((bn[mx]+bn[mx+1])/2,2))
                    res = get_median_with_error(cnt,bn)
                    ax_ij.set_yticks([])
                    if labels is not None:
                        if udm is None:
                            udm = ["\"" for i in range(len(labels))]
                        elif type(udm)==str:
                            udm = [udm for i in range(len(labels))]
                        ax_ij.set_title("$"+labels[i].replace("$","")+"$" +" = "+ res+" "+udm[i])
                    
                else:
                    if [i,j]==[1,0]:
                        mat = prob3D.sum(axis=2)#/np.sum(Combined_PDF.sum(axis=2))
                        if labels is not None:
                            ax_ij.set_ylabel(labels[1])
                    elif [i,j]==[2,0]:
                        mat = prob3D.sum(axis=1)#/np.sum(Combined_PDF.sum(axis=1))
                        if labels is not None:
                            ax_ij.set_xlabel(labels[0])
                            ax_ij.set_ylabel(labels[2])
                    elif [i,j]==[2,1]:
                        mat = prob3D.sum(axis=0)#/np.sum(Combined_PDF.sum(axis=0) )
                        ax_ij.set_xlabel(labels[1])
                    ax_ij.contour(bin_coord[j],bin_coord[i],mat.T,levels=contour_levels,alpha=alpha)#**contourf_kwargs) #or ax_ij.contourf(bin_coord[j],bin_coord[i],mat.T)
                    #ax_ij.pcolormesh(bin_coord[j],bin_coord[i],mat.T,shading='gouraud',)
                    ax_ij.set_xticks(ticks[j])
                    ax_ij.set_yticks(ticks[i])
            else:
                ax_ij.axis("off")
    
    plt.tight_layout()

    return plt


# In[ ]:


def plot_probability3D_KDE(KDE,Positions,labels=None,udm=None,figsize=(12,12),contour_levels=30,colors=base_colors,alpha=1):
    fg,ax= plt.subplots(3,3,figsize=figsize) 
    #bin_coord = [[(bins[j][i]+bins[j][i+1])/2. for i in range(0,len(bins[j])-1)]for j in range(len(bins))]
    #bin_densities = marginalise_prob(prob3D,bins) 
    ticks = []
    """for i in range(len(bins)):
        rng = bins[i][-1]- bins[i][0]
        stp = bins[i][1] - bins[i][0]
        ticks_i=[]
        for j in range(len(bins[i])):
            ticks_i.append(np.round(bins[i][0]+j*stp,1))
        ticks.append(ticks_i)"""
        
    for i in range(3):
        for j in range(3):
            ax_ij = ax[i][j]
            if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                if i==j:
                    if labels is not None:
                        if i==0:
                            ax_ij.set_ylabel(labels[i])
                        elif i==2:
                            ax_ij.set_xlabel(labels[i])

                            
                    x = Positions[i][::len(KDE)]
                    y = KDE[i]/np.sum(KDE[i])                            
                    ax_ij.plot(x,y)

                    #ax_ij.set_xticks(ticks[j])
                    #ax_ij.set_ylim(0,max(cnt)+.2*max(cnt))
                    mx = np.where(y==max(y))[0][0]
                    #res = str(np.round((bn[mx]+bn[mx+1])/2,2))
                    #res = get_median_with_error(cnt,bn)
                    ax_ij.set_yticks([])
                    if labels is not None:
                        if udm is None:
                            udm = ["\"" for i in range(len(labels))]
                        elif type(udm)==str:
                            udm = [udm for i in range(len(labels))]
                        ax_ij.set_title("$"+labels[i].replace("$","")+udm[i])#+"$" +" = "+ res+" "+udm[i])
                    
                else:
                    if [i,j]==[1,0]:
                        #mat = prob3D.sum(axis=2)#/np.sum(Combined_PDF.sum(axis=2))
                        mat = np.array(KDE[:3])
                        if labels is not None:
                            ax_ij.set_ylabel(labels[1])
                    elif [i,j]==[2,0]:
                        mat = np.array([KDE[0],KDE[2]])

                        if labels is not None:
                            ax_ij.set_xlabel(labels[0])
                            ax_ij.set_ylabel(labels[2])
                    elif [i,j]==[2,1]:
                        mat = np.array(KDE[1:])
                        ax_ij.set_xlabel(labels[1])
                    #TEST
                    xx = Positions[i][::len(KDE)]
                    yy = Positions[j][::len(KDE)] 
                    ax_ij.contour(xx, yy, mat, levels=4,alpha=alpha,colors=colors[i])

                    #ax_ij.pcolormesh(bin_coord[j],bin_coord[i],mat.T,shading='gouraud',)
                    #ax_ij.set_xticks(ticks[j])
                    #ax_ij.set_yticks(ticks[i])
            else:
                ax_ij.axis("off")
    
    plt.tight_layout()

    return plt


# In[ ]:


def plot_model_WS(modelPlot,savefig_path,v_min,v_max,res_min,res_max):
    
    f, axes = plt.subplots(1, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.data_plot(ax=axes[0], v_min=v_min, v_max=v_max)
    modelPlot.model_plot(ax=axes[1], v_min=v_min, v_max=v_max)
    modelPlot.normalized_residual_plot(ax=axes[2],v_min=res_min, v_max=res_max)
    plt.savefig(savefig_path+"/results_Model.png")
    f, axes = plt.subplots(1, 2, figsize=(10, 8), sharex=False, sharey=False)
    modelPlot.convergence_plot(ax=axes[0], v_max=1)
    modelPlot.magnification_plot(ax=axes[ 1])
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.05)
    plt.savefig(savefig_path+"/results_MassModel.png")

    f, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,1], text='All components', source_add=True, lens_light_add=True,
                                 unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,1], text='All components convolved', source_add=True, \
                                 lens_light_add=True, point_source_add=True, v_min=v_min, v_max=v_max)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig(savefig_path+"/results_LightComponents.png")
    plt.close()

    
def plot_model(modelPlot,savefig_path,v_min,v_max,res_min,res_max):
    
    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.data_plot(ax=axes[0,0], v_min=v_min, v_max=v_max)
    modelPlot.model_plot(ax=axes[0,1], v_min=v_min, v_max=v_max)
    modelPlot.normalized_residual_plot(ax=axes[0,2],v_min=res_min, v_max=res_max)
    modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=1000)
    modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
    modelPlot.magnification_plot(ax=axes[1, 2])
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig(savefig_path+"/results_Model.png")
    plt.close()

    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True,
                                 unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, \
                                 lens_light_add=True, point_source_add=True, v_min=v_min, v_max=v_max)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig(savefig_path+"/results_LightComponents.png")
    plt.close()


# In[ ]:


from statistical_tools import estimate_sigma, estimate_median
def get_median_with_error(prob,bins):
    median = estimate_median([prob],[bins])[0]
    sigmas = estimate_sigma([prob],[bins],median=[median],averaged=False)[0]
    rounding = int(1-np.log10(max(sigmas)))
    return str(np.round(median,rounding))+"$_{"+str(np.round(sigmas[0],rounding))+"}^{"+str(np.round(sigmas[1],rounding))+"}$"


# In[ ]:


import scipy.stats as st
# Discontinued
"""
#from corner import quantile

def my_corner(samples,samples_names=None,labels=None,udm=None,figsize=(12,12),contour_levels=4,colors=base_colors,alpha=.3):
    
    if samples_names is None or len(samples_names)<len(samples):
        samples_names = ["" for i in samples]
        
    fg,ax= plt.subplots(3,3,figsize=figsize) 
    
    # bin determination as in Combined_PDF
    n_bins = []
    ticks  = []
    for i in range(len(samples[0])):
        sT    = np.array(samples,dtype=object)[:,i] 
        top_i = min([len(k) for k in sT])
        n_bins_i = 3*int(round( 1 + 3.322*np.log10(top_i)))+2
        n_bins.append(n_bins_i)
        
        max_data_i = max([max(k) for k in sT])
        min_data_i = min([min(k) for k in sT])
        diff_min_max = (max_data_i-min_data_i)
        dbin_i = diff_min_max/n_bins_i
        ticks_i=[np.round(min_data_i,2)]
        for nb in range(n_bins_i-1):
            ticks_i.append(np.round(ticks_i[-1]+dbin_i,2))
        ticks.append(ticks_i)
        
        
    for j,smpl in enumerate(samples):
        col= colors[j]

        MCMC  = np.array(smpl)
        print("Sample ",samples_names[j])
        for i in range(3):
            for j in range(3):
                ax_ij = ax[i][j]
                if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                    if i==j:
                        if labels is not None:
                            if i==0:
                                ax_ij.set_ylabel(labels[i])
                            elif i==2:
                                ax_ij.set_xlabel(labels[i])
                        #if type(bins) is int:
                        #    bins=[bins for j in range(3)]
                        ax_ij.hist(MCMC[i],bins=n_bins[i],color=col,alpha=alpha,density=True)
                        vmin,res,vmax = quantile(MCMC[i],q=[0.16, 0.5, 0.84])
                        err_min = res-vmin
                        err_max = vmax-res
                        rndto = 0
                        while np.round(min([err_min,err_max]),rndto)<=0.:
                            rndto+=1
                        res_str = "$"+str(np.round(res,rndto))+"^{+"+str(np.round(err_min,rndto))+"}_{-"+str(np.round(err_max,rndto))+"} \ "+udm[i]+"$"
                        ax_ij.set_xticks(ticks[j])
                        ax_ij.set_yticks([])
                        if labels is not None:
                            if udm is None:
                                udm = ["\"" for i in range(len(labels))]
                            elif type(udm)==str:
                                udm = [udm for i in range(len(labels))]
                            ax_ij.set_title("$"+labels[i].replace("$","")+"$" +" = "+ res_str)
                    else:
                        if [i,j]==[1,0]:
                            mat = MCMC[:3] 
                            if labels is not None:
                                ax_ij.set_ylabel(labels[1])
                        elif [i,j]==[2,0]:
                            mat = np.array([MCMC[0],MCMC[2]])
                            if labels is not None:
                                ax_ij.set_xlabel(labels[0])
                                ax_ij.set_ylabel(labels[2])
                        elif [i,j]==[2,1]:
                            mat = MCMC[1:]
                            ax_ij.set_xlabel(labels[1])
                        xmin,xmax = min(mat[0]),max(mat[0])
                        ymin,ymax = min(mat[1]),max(mat[1])
                        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                        positions = np.vstack([xx.ravel(), yy.ravel()])
                        #mat_small =int(len(mat[0])/2)
                        n_points = 100000
                        #values = np.vstack([mat[0][mat_small:mat_small+n_points] , mat[1][mat_small:mat_small+n_points]])
                        values = np.vstack([mat[0][-n_points:] , mat[1][-n_points:]])
                        kernel = st.gaussian_kde(values)
                        f = np.reshape(kernel(positions).T, xx.shape)
                        ax_ij.contour(xx, yy, f, levels=contour_levels,alpha=alpha,colors=col)
                        ax_ij.set_xticks(ticks[j])
                        ax_ij.set_yticks(ticks[i])
                else:
                    if [i,j]!=[2,0]:
                        ax_ij.axis("off")
    axdel=ax[0][2]
    for i,set_name in enumerate(samples_names):
        axdel.scatter([1],[1],label=set_name.replace("_"," "),color=colors[i])
    axdel.scatter([1],[1], color="w")
    axdel.legend()
    axdel.axis("off")


    plt.tight_layout()
    return plt
"""


# In[ ]:


from tools import *

def my_corner_general(samples,samples_names=None,labels=None,udm=None,figsize=None,contour_levels=4,
                      colors=base_colors,alpha=.3):
    if figsize is None:
        figsize= (3*len(samples[0]),3*len(samples[0]))
    if samples_names is None or len(samples_names)<len(samples):
        samples_names = ["" for i in samples]
        
    fg,ax= plt.subplots(len(samples[0]),len(samples[0]),figsize=figsize) 
    
    # bin determination as in Combined_PDF
    n_bins = []
    ticks  = []
    for i in range(len(samples[0])):
        sT    = np.array(samples,dtype=object)[:,i] 

        top_i = min([len(k) for k in sT])
        n_bins_i = 3*int(round( 1 + 3.322*np.log10(top_i)))+2
        
        n_bins.append(n_bins_i)
        
        max_data_i = max([max(k) for k in sT])
        min_data_i = min([min(k) for k in sT])
        diff_min_max = (max_data_i-min_data_i)
        n_ticks = 6
        step    = diff_min_max/n_ticks
        
        rnd_ticks = 0
        while True:
            if np.round(step,rnd_ticks)>0:
                break
            else:
                rnd_ticks+=1 
        if rnd_ticks<=2:
            rnd_ticks += 1
        
        ticks_i = [np.round(min_data_i,rnd_ticks)]
        for i in range(1,n_ticks):
            ticks_i.append(np.round(ticks_i[i-1]+step,rnd_ticks))
            
        ticks.append(ticks_i)

        
        
    for k,smpl in enumerate(samples): # stg
        col= colors[k]

        print("Sample ",samples_names[k])
        MCMC  = np.array(smpl) 
        
        for i in range(len(samples[0])): # prm
            for j in range(len(samples[0])):
                ax_ij = ax[i][j]
                if j<=i:
                    if i==j:
                        ax_ij.hist(MCMC[i],bins=n_bins[i],color=col,alpha=alpha,density=True)
                        #if k==0:
                        ax_ij.set_xticks(ticks[i])
                        ax_ij.set_yticks([])
                        if labels is not None:
                            ax_ij.set_title(labels[i])
                            if i==len(samples[0])-1 :
                                ax_ij.set_xlabel(labels[i])
                    else:
                        mat = [MCMC[j],MCMC[i]]
                        if labels is not None and k==0:
                            if j==0:
                                ax_ij.set_ylabel(labels[i])
                            if i==len(samples[0])-1:
                                ax_ij.set_xlabel(labels[j])
                        xmin,xmax = min(mat[0]),max(mat[0])
                        ymin,ymax = min(mat[1]),max(mat[1])
                        xx, yy    = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
                        positions = np.vstack([xx.ravel(), yy.ravel()])
                        #mat_small =int(len(mat[0])/2)
                        n_points = 100000
                        #values = np.vstack([mat[0][mat_small:mat_small+n_points] , mat[1][mat_small:mat_small+n_points]])
                        values = np.vstack([mat[0][-n_points:] , mat[1][-n_points:]])
                        kernel = st.gaussian_kde(values)
                        f = np.reshape(kernel(positions).T, xx.shape)
                        ax_ij.contour(xx, yy, f, levels=contour_levels,alpha=alpha,colors=col)
                        #if k==0:
                        ax_ij.set_yticks(ticks[i])
                        ax_ij.set_xticks(ticks[j])
                        if j!=0:
                            ax_ij.set_yticklabels([" " for n in range(len(ticks[i]))])
                    if i!=len(samples[0])-1 and k==0:                    
                        ax_ij.set_xticklabels([" " for n in range(len(ticks[j]))])
                else:
                    ax_ij.axis("off")
    axdel=ax[0][-1]
    for i,set_name in enumerate(samples_names):
        axdel.scatter([1],[1],label=strip_setting_name(set_name).replace("_"," "),color=colors[i])
    axdel.scatter([1],[1], color="w")
    axdel.legend()
    axdel.axis("off")
    plt.tight_layout()
    plt.subplots_adjust(wspace=.06, hspace=.06)

    return plt


# In[ ]:


def overplot_probability3D(prob3D_list,bins_list,labels=None,labels_list=None,udm=None,figsize=(12,12),contour_levels=30,colors=[["r","darkorange"],["b","royalblue"],["g","lime"]]):
    fg,ax= plt.subplots(3,3,figsize=figsize) 
    alpha=0.3
    list_index = 0
    for prob3D,bins in zip(prob3D_list,bins_list):
        bin_coord = [[(bins[j][i]+bins[j][i+1])/2. for i in range(0,len(bins[j])-1)]for j in range(len(bins))]
        bin_densities = marginalise_prob(prob3D,bins) 
        ticks = []
        for i in range(len(bins)):
            rng = bins[i][-1]- bins[i][0]
            stp = bins[i][1] - bins[i][0]
            ticks_i=[]
            for j in range(len(bins[i])):
                ticks_i.append(np.round(bins[i][0]+j*stp,1))
            ticks.append(ticks_i)

        for i in range(3):
            for j in range(3):
                ax_ij = ax[i][j]
                if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                    if i==j:
                        if labels is not None:
                            if i==0:
                                ax_ij.set_ylabel(labels[i])
                            elif i==2:
                                ax_ij.set_xlabel(labels[i])
                        cnt = bin_densities[i]
                        bn  = bins[i]
                        for k in range(len(cnt)):
                            lbl=None
                            if k==len(cnt)-1:
                                lbl = labels_list[list_index]
                            ax_ij.fill_between([bn[k],bn[k+1]],[0],cnt[k],color=colors[i][list_index],alpha=alpha,label=lbl)
                        ax_ij.set_xticks(ticks[j])
                        ax_ij.set_yticks([])
                        ax_ij.set_ylim(0,max(cnt)+.2*max(cnt))
                        if labels is not None:
                            if udm is None:
                                udm = ["\"" for i in range(len(labels))]
                            elif type(udm)==str:
                                udm = [udm for i in range(len(labels))]
                            ax_ij.set_title("$"+labels[i].replace("$","")+"$")
                        ax_ij.legend()
                    else:
                        if [i,j]==[1,0]:
                            mat = prob3D.sum(axis=2)#/np.sum(Combined_PDF.sum(axis=2))
                            if labels is not None:
                                ax_ij.set_ylabel(labels[1])
                        elif [i,j]==[2,0]:
                            mat = prob3D.sum(axis=1)#/np.sum(Combined_PDF.sum(axis=1))
                            if labels is not None:
                                ax_ij.set_xlabel(labels[0])
                                ax_ij.set_ylabel(labels[2])
                        elif [i,j]==[2,1]:
                            mat = prob3D.sum(axis=0)#/np.sum(Combined_PDF.sum(axis=0) )
                            ax_ij.set_xlabel(labels[1])
                        ax_ij.contour(bin_coord[j],bin_coord[i],mat.T,levels=contour_levels,alpha=alpha)#**contourf_kwargs) #or ax_ij.contourf(bin_coord[j],bin_coord[i],mat.T)
                        #ax_ij.pcolormesh(bin_coord[j],bin_coord[i],mat.T,shading='gouraud',)
                        ax_ij.set_xticks(ticks[j])
                        ax_ij.set_yticks(ticks[i])
                else:
                    ax_ij.axis("off")
        list_index+=1

    plt.tight_layout()

    return plt

