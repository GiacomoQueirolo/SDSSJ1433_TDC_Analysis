#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
import matplotlib.pyplot as plt

from Utils.tools import *
from Utils.statistical_tools import marginalise_prob
base_colors = ["r","b","g","y","k","m","c","darkorange","darkviolet","lawngreen","violet"] 

def plot_probability3D(prob3D,bins,labels=None,udm=None,figsize=(16,16),num_contour_levels=6,colors=base_colors,alpha=1,title="",fnt = 20,cont_sig=True):
    plt.rcParams['xtick.labelsize'] = fnt
    plt.rcParams['ytick.labelsize'] = fnt 
    plt.rcParams['font.size'] = fnt
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('font', size=fnt)          # controls default text sizes
    plt.rc('axes', titlesize=fnt)     # fontsize of the axes title
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fnt)    # legend fontsize
    
    fg,ax= plt.subplots(3,3,figsize=figsize) 
    bin_coord = [[(bins[j][i]+bins[j][i+1])/2. for i in range(0,len(bins[j])-1)]for j in range(len(bins))]
    bin_densities = marginalise_prob(prob3D,bins) 
    ticks = []
    for i in range(len(bins)):
        rng = bins[i][-1]- bins[i][0]
        stp = bins[i][1] - bins[i][0]
        #for j in range(0,len(bins[i]),):
        #    ticks_i.append(np.round(bins[i][0]+j*stp,1))
        #ticks_i =np.arange(min(ticks_i),max(ticks_i),.1)
        n_ticks = 4
        ln_stp  = rng/n_ticks 
        int_rnd = int(abs(np.floor(np.log10(abs(ln_stp)))))
        ticks_i=[]
        for j in range(n_ticks+1):
            ticks_i.append(np.round(bins[i][0]+j*(ln_stp),int_rnd+1))
        ticks.append(ticks_i) 
        
    for i in range(3):
        for j in range(3):
            ax_ij = ax[i][j]
            if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                if i==j:
                    if labels is not None:
                        if i==2:
                            ax_ij.set_xlabel(labels[i])
                    cnt = bin_densities[i]
                    bn  = bins[i]
                    for k in range(len(cnt)):
                        ax_ij.fill_between([bn[k],bn[k+1]],[0],cnt[k],color=colors[i],alpha=alpha)

                    ax_ij.set_ylim(0,max(cnt)+.2*max(cnt))
                    res = get_median_with_error(cnt,bn)
                    ax_ij.set_xticks(ticks[j])
                    if j!=2:
                        ax_ij.set_xticklabels([""]*len(ticks[j]))
                    else:
                        ax_ij.set_xticklabels(ticks[j],rotation=45)
                    ax_ij.set_yticks([])
                    if labels is not None:
                        if udm is None:
                            udm = [r"$arcsec^2$" for i in range(len(labels))]
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
                        if labels is not None:
                            ax_ij.set_xlabel(labels[1])
                    if not cont_sig:
                        contour_levels = [((i+1)*mat.max()/float(num_contour_levels))  for i in range(num_contour_levels)]
                        #lab_cnt = None
                    else:
                        contour_levels = [(100-n)*mat.max()/100  for n in [99.7,95,68] ]
                        #lab_cnt = r"1,2 and 3 $\sigma$ confidence regions"
                    ax_ij.contour(bin_coord[j],bin_coord[i],mat.T,levels=contour_levels,alpha=alpha)
                    ax_ij.set_xticks(ticks[j]) 
                    ax_ij.set_yticks(ticks[i])  
                    if [i,j]==[1,0]:
                        ax_ij.set_xticklabels([""]*len(ticks[j]))
                        ax_ij.set_yticklabels(ticks[i],rotation=45)
                    elif [i,j]==[2,1]:
                        ax_ij.set_yticklabels([""]*len(ticks[i]))
                        ax_ij.set_xticklabels(ticks[j],rotation=45)
                    else:
                        ax_ij.set_yticklabels(ticks[i],rotation=45)
                        ax_ij.set_xticklabels(ticks[j],rotation=45)
                    

   
            else:
                ax_ij.axis("off")
    if title:
        fg.suptitle(title)
    plt.tight_layout(pad=1.3)
    plt.subplots_adjust(wspace=0.1,
                    hspace=0.1)
    return plt



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
                    
                    #mx = np.where(y==max(y))[0][0]
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
 
def plot_model_WS(modelPlot,savefig_path,v_min,v_max,res_min,res_max,band=(0,""),verbose=True):
    band_index,filter = band
    if filter!="":
        filter="_"+filter
    f, axes = plt.subplots(1, 3, figsize=(16, 8), sharex=False, sharey=False)
    modelPlot.data_plot(ax=axes[0], band_index=band_index,v_min=v_min, v_max=v_max)
    modelPlot.model_plot(ax=axes[1], band_index=band_index, v_min=v_min, v_max=v_max)
    modelPlot.normalized_residual_plot(ax=axes[2], band_index=band_index,v_min=res_min, v_max=res_max)
    savename= f"{savefig_path}/results_Model{filter}.png"
    if verbose:
        print(f"Saved model plot {savename}") 
    plt.savefig(savename)
    
    f, axes = plt.subplots(1, 2, figsize=(10, 8), sharex=False, sharey=False)
    modelPlot.convergence_plot(ax=axes[0], band_index=band_index, v_max=1)
    modelPlot.magnification_plot(ax=axes[ 1], band_index=band_index)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.05)
    savename = f"{savefig_path}/results_MassModel{filter}.png"
    if verbose:
        print(f"Saved mass model plot {savename}") 
    plt.savefig(savename)
    plt.close()

    f, axes = plt.subplots(2, 2, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.decomposition_plot(ax=axes[0,0], band_index=band_index, text='Lens light', lens_light_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,0], band_index=band_index, text='Lens light convolved', lens_light_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,1], band_index=band_index, text='All components', source_add=True, lens_light_add=True,
                                 unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,1], band_index=band_index, text='All components convolved', source_add=True, \
                                 lens_light_add=True, point_source_add=True, v_min=v_min, v_max=v_max)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)

    savename = f"{savefig_path}/results_LightComponents{filter}.png"
    if verbose:
        print(f"Saved light model plot {savename}") 
    plt.savefig(savename)
    plt.close()

    
def plot_model(modelPlot,savefig_path,v_min,v_max,res_min,res_max,band=(0,""),verbose=True):
    band_index,filter = band
    if filter!="":
        filter="_"+filter
    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.data_plot(ax=axes[0,0],band_index=band_index, v_min=v_min, v_max=v_max)
    modelPlot.model_plot(ax=axes[0,1],band_index=band_index, v_min=v_min, v_max=v_max)
    modelPlot.normalized_residual_plot(ax=axes[0,2],band_index=band_index,v_min=res_min, v_max=res_max)
    modelPlot.source_plot(ax=axes[1, 0],band_index=band_index, deltaPix_source=0.01, numPix=1000)
    modelPlot.convergence_plot(ax=axes[1, 1],band_index=band_index, v_max=1)
    modelPlot.magnification_plot(ax=axes[1, 2],band_index=band_index)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    savename = f"{savefig_path}/results_MassModel{filter}.png"
    if verbose:
        print(f"Saved mass model plot {savename}") 
    plt.savefig(savename)
    plt.close()

    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    modelPlot.decomposition_plot(ax=axes[0,0],band_index=band_index, text='Lens light', lens_light_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,0],band_index=band_index, text='Lens light convolved', lens_light_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,1],band_index=band_index, text='Source light', source_add=True, unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,1],band_index=band_index, text='Source light convolved', source_add=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[0,2],band_index=band_index, text='All components', source_add=True, lens_light_add=True,
                                 unconvolved=True, v_min=v_min, v_max=v_max)
    modelPlot.decomposition_plot(ax=axes[1,2],band_index=band_index, text='All components convolved', source_add=True, \
                                 lens_light_add=True, point_source_add=True, v_min=v_min, v_max=v_max)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    savename = f"{savefig_path}/results_LightComponents{filter}.png"
    if verbose:
        print(f"Saved light model plot {savename}") 
    plt.savefig(savename)
    plt.close()



from Utils.statistical_tools import estimate_sigma, estimate_median
def get_median_with_error(prob,bins,ret_str=True):
    median = estimate_median([prob],[bins])[0]
    sigmas = estimate_sigma([prob],[bins],median=[median],averaged=False)[0]
    rounding = int(1-np.log10(max(sigmas)))
    if ret_str:
        sig_up_rnd,sig_down_rnd = str(np.round(sigmas[0],rounding)),str(np.round(sigmas[1],rounding))
        if sig_up_rnd==sig_down_rnd:
            return str(np.round(median,rounding))+"$\pm$"+sig_up_rnd
        return str(np.round(median,rounding))+"$_{"+sig_down_rnd+"}^{"+sig_up_rnd+"}$"
    else:
        return median,sigmas

def _print_3Dmarg_res(Combined_PDF,Combined_bins):
    bin_densities = marginalise_prob(Combined_PDF,Combined_bins) 
    for i in range(3):
        bn  =  Combined_bins[i]
        cnt = bin_densities[i]
        res = get_median_with_error(cnt,bn)
        print(res)
    return 0

import scipy.stats as st
""" # Discontinued
from corner import quantile


#def my_corner(samples,samples_names=None,labels=None,udm=None,figsize=(12,12),contour_levels=4,colors=base_colors,alpha=.3):
    
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
        ticks_i=[np.round(min_data_i,1)]
        for nb in range(n_bins_i-1):
            ticks_i.append(np.round(ticks_i[-1]+dbin_i,1))
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
                        ax_ij.hist(MCMC[i],bins=n_bins[i],color=col,alpha=alpha,density=True,stacked=True)
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
                            ax_ij.set_title(labels[i])

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

#samples shape = setting, param, steps

def my_corner_general(samples,samples_names=None,labels=None,udm=None,figsize=None,contour_levels=4,
                      colors=base_colors,alpha=.3):
    fnt = 30
    plt.rcParams['xtick.labelsize'] = fnt
    plt.rcParams['ytick.labelsize'] = fnt 
    plt.rcParams['font.size'] = fnt
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('font', size=fnt)          # controls default text sizes
    plt.rc('axes', titlesize=fnt)     # fontsize of the axes title
    plt.rc('axes', labelsize=fnt)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=fnt)    # fontsize of the tick labels
    plt.rc('legend', fontsize=fnt)    # legend fontsize
    
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
                        ax_ij.hist(MCMC[i],bins=n_bins[i],color=col,alpha=alpha,density=True,stacked=True)
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

def averaged_plot(ax,prm_steps,col="b", num_average=100,plot_scatter=False,param_name=None,renorm=False):
    # general plot for the behaviour of non-linear solver for 1 param
    """
    :param prm_steps: parameters sampled 2d numpy array
    :param num_average: number of samples to average (should coincide with the number of samples in the emcee process)
    :return:
    """
    prm_steps = np.array(prm_steps)
    num_average = int(num_average)
    num_steps   = len(prm_steps)
    n_points = int((num_steps - num_steps % num_average) / num_average)
    sliced_steps = prm_steps[:int(n_points * num_average)].reshape(n_points, num_average)
    steps_averaged = np.average(sliced_steps, axis=1)
    if renorm:
        end_point = np.mean(steps_averaged)
        steps_renormed = (steps_averaged - end_point) / np.std(steps_averaged)
    else:
        steps_renormed = steps_averaged
    x=np.arange(0,len(steps_renormed)) 
    ax.plot(x,steps_renormed, label=param_name,color=col)
    mn = min(steps_renormed)
    fac_mn = 1-0.002
    if mn<0:
        fac_mn = 1+0.002
    mx = max(steps_renormed)
    fac_mx = 1+0.002
    if mx<0:
        fac_mx = 1-0.002
    ax.set_ylim(mn*fac_mn,mx*fac_mx)

    if plot_scatter:
        # MOD_SCATTER
        sigma_up,sigma_down = [],[]
        for j in range(len(steps_averaged)):
            sg = np.std(sliced_steps[j])
            sigma_up.append(steps_averaged[j]+sg)
            sigma_down.append(steps_averaged[j]-sg)
        if renorm:
            sigma_up_ren   = (sigma_up-end_point)/np.std(steps_averaged)
            sigma_down_ren = (sigma_down-end_point)/np.std(steps_averaged)
        else:
            sigma_up_ren   = sigma_up
            sigma_down_ren = sigma_down
        ax.fill_between(x,sigma_up_ren,sigma_down_ren,alpha=0.3,color=col,label=r"1-$\sigma$ scatter")            
    ax.legend()
    return ax
"""
from Data.input_data import *
from Utils.get_res import get_kwres 
from lenstronomy.Plots.model_plot import ModelPlot


def get_ModelPlot(setting,kwargs_data=None, kwargs_psf=None, kwargs_numerics=None,kwargs_model=None,\
                  kwargs_results=None,image_likelihood_mask_list=None,arrow_size=0.02, cmap_string="gist_heat",withmask=True):
    if kwargs_data is None:
        kwargs_data,mask = init_kwrg_data(setting,saveplots=False,return_mask=True)
    if kwargs_psf is None:
        kwargs_psf = init_kwrg_psf(setting,saveplots=False)
    if kwargs_numerics is None:
        kwargs_numerics = init_kwrg_numerics(setting)
    if kwargs_model is None:
        kwargs_model = get_kwargs_model(setting)
    if kwargs_results is None:
        kwargs_results = get_kwres(setting)["kwargs_results"]
    if image_likelihood_mask_list is None and withmask:
        try:
            image_likelihood_mask_list = [np.array(mask).tolist()]
        except:
            _,mask  = init_kwrg_data(setting,saveplots=False,return_mask=True)
            image_likelihood_mask_list = [np.array(mask).tolist()]
    else:
        image_likelihood_mask_list = []
    multi_band = [[kwargs_data,kwargs_psf,kwargs_numerics]]
    modelplot = ModelPlot(multi_band,kwargs_model,kwargs_results,\
                image_likelihood_mask_list=image_likelihood_mask_list,arrow_size=arrow_size,cmap_string=cmap_string)
    return modelplot
"""
def _get_ModelPlot(kwargs_data, kwargs_psf, kwargs_numerics,kwargs_model,\
                  kwargs_results,image_likelihood_mask_list,arrow_size=0.02, cmap_string="gist_heat"):
    multi_band = [[kwargs_data,kwargs_psf,kwargs_numerics]]
    modelplot = ModelPlot(multi_band,kwargs_model,kwargs_results,\
                image_likelihood_mask_list=image_likelihood_mask_list,arrow_size=arrow_size,cmap_string=cmap_string)
    return modelplot
