#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#With argparse which should help make it easier


# In[1]:


import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

from Utils.tools import *
from Utils.get_res import *
from Plots.plotting_tools import base_colors,averaged_plot

def my_plot_mcmc_behaviour(ax, samples_mcmc, param_mcmc, col="b", num_average=100):
    """
    plots the MCMC behaviour and looks for convergence of the chain
    :param samples_mcmc: parameters sampled 2d numpy array
    :param param_mcmc: list of parameters
    :param num_average: number of samples to average (should coincide with the number of samples in the emcee process)
    :return:
    """
    num_samples = len(samples_mcmc[:, 0])
    num_average = int(num_average)
    n_points = int((num_samples - num_samples % num_average) / num_average)
    for i, param_name in enumerate(param_mcmc):
        samples = samples_mcmc[:, i]
        sliced_sample = samples[:int(n_points * num_average)].reshape(n_points, num_average)
        samples_averaged = np.average(sliced_sample, axis=1)
        end_point = np.mean(samples_averaged)
        samples_renormed = (samples_averaged - end_point) / np.std(samples_averaged)
        # MOD_SCATTER
        sigma_up,sigma_down = [],[]
        for j in range(len(samples_averaged)):
            sg = np.std(sliced_sample[j])
            sigma_up.append(samples_averaged[j]+sg)
            sigma_down.append(samples_averaged[j]-sg)
        sigma_up_ren   = (sigma_up-end_point)/np.std(samples_averaged)
        sigma_down_ren = (sigma_down-end_point)/np.std(samples_averaged)
        x=np.arange(0,len(samples_renormed))
        ax.fill_between(x,sigma_up_ren,sigma_down_ren,alpha=0.3,color=col,label=r"1-$\sigma$ scatter")
        ax.plot(x,samples_renormed, label=param_name,color=col)
        mn = min(samples_renormed)
        fac_mn = 1-0.002
        if mn<0:
            fac_mn = 1+0.002
        mx = max(samples_renormed)
        fac_mx = 1+0.002
        if mx<0:
            fac_mx = 1-0.002
        ax.set_ylim(mn*fac_mn,mx*fac_mx)
    ax.legend()
    return ax

def my_plot_mcmc_simil_behaviour(ax, samples_mcmc, param_mcmc,col="b", num_average=100,boundary=None):
    """
    plots the MCMC behaviour and looks for convergence of the chain
    :param samples_mcmc: parameters sampled 2d numpy array
    :param param_mcmc: list of parameters
    :param num_average: number of samples to average 
        (should coincide with the number of samples in the emcee process)
    :return:
    """
    num_samples = len(samples_mcmc[:, 0])
    num_average = int(num_average)
    n_points = int((num_samples - num_samples % num_average) / num_average)
    for i, param_name in enumerate(param_mcmc):
        samples = samples_mcmc[:, i]
        sliced_sample = samples[:int(n_points * num_average)].reshape(n_points, num_average)
        samples_averaged = np.average(sliced_sample, axis=1)
        # MOD_SCATTER
        sigma_up,sigma_down = [],[] 
        for j in range(len(samples_averaged)):
            sg = np.std(sliced_sample[j])
            sigma_up.append(samples_averaged[j]+sg)
            sigma_down.append(samples_averaged[j]-sg)
        x=np.arange(0,len(samples_averaged))
        
        ax.fill_between(x,sigma_up,sigma_down,alpha=0.3,color=col)
        ax.plot(x,samples_averaged,label=param_name,color=col)
        mn = min(samples_averaged)
        fac_mn = 1-0.002
        if mn<0:
            fac_mn = 1+0.002
        mx = max(samples_averaged)
        fac_mx = 1+0.002
        if mx<0:
            fac_mx = 1-0.002
        ax.set_ylim(mn*fac_mn,mx*fac_mx)

        #MOD_BOUNDARIES
        if boundary != None and "ra_" not in param_mcmc[0] and "dec_" not in param_mcmc[0]: 
            lim_up = max(samples_averaged)
            lim_dwn = min(samples_averaged)
            diff_lim = lim_up-lim_dwn
            lim_up +=diff_lim/14.
            lim_dwn -=diff_lim/14.
            if min(boundary)>=(lim_dwn-4*diff_lim/14.):
                ax.hlines(boundary[0],0,len(samples_averaged),colors="r",linestyles='dashed',label="lower bound")
            elif max(boundary)<=(lim_up+4*diff_lim/14.):
                ax.hlines(boundary[1],0,len(samples_averaged),colors="r",linestyles='dashed',label="upper bound")
                
    ax.legend()
    return ax

def find_bound(param_name, setting,num=None):
    setting  = get_setting_module(setting).setting()
    if "lens_light" in param_name:
        lens_light_params = setting.lens_light_params 
        n_lens = int(param_name[-1])
        lower  = lens_light_params[3][n_lens]       
        upper  = lens_light_params[4][n_lens]
        prm    = param_name.replace("_lens_light"+str(n_lens),"")
    elif "lens" in param_name:
        lens_params = setting.lens_params 
        n_lens = int(param_name[-1])
        lower  = lens_params[3][n_lens]       
        upper  = lens_params[4][n_lens]
        prm    = param_name.replace("_lens"+str(n_lens),"")
    elif "source_light0" in param_name:
        source_params = setting.source_params
        lower = source_params[3][0]       
        upper = source_params[4][0]        
        prm   = param_name.replace("_source_light0","")
    elif "image" in param_name:
        if "ra_" in param_name:
            radec="ra_image"
        else:
            radec="dec_image"
        ps_params = setting.ps_params
        lower     = ps_params[3][0][radec]
        upper     = ps_params[4][0][radec]
        prm       = num
    else:
        raise ValueError("Unrecognised parameter name:"+param_name)
    return lower[prm],upper[prm]


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the MCMC behaviour plots of ALL parameter separately for the given filter model")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-al", "--all_logL", dest="all_logL", default=True, action="store_false",
                        help="Plot every point of the log likelihood and its final histogram")
    parser.add_argument('SETTING_FILE', help="setting file to consider")

    args     = parser.parse_args()
    cut_mcmc = args.cut_mcmc
    all_logL = args.all_logL

    present_program(sys.argv[0])

    setting      = get_setting_module(args.SETTING_FILE).setting()
    backup_path  = "backup_results/" 
    savefig_path = get_savefigpath(setting,backup_path)+"/MCMC_bhv/"
    mkdir(savefig_path)


    param_mcmc   = get_mcmc_prm(setting,backup_path)
    samples_mcmc = get_mcmc_smpl(setting,backup_path)[cut_mcmc:]
    mcmc_fermat  = get_mcmc_fermat(setting,backup_path)[cut_mcmc:]
    try:
        mcmc_logL = get_mcmc_logL(setting,backup_path)
        logL_go =True
    except FileNotFoundError:
        print(mcmc_logL_file+" not found")
        logL_go = False

    mkdir(savefig_path+"/simil_mcmc_bhv")
        
    colors = base_colors
    j=0
    n_ra=0
    n_dec=0
    for i in range(len(param_mcmc)):  
        f, axes = plt.subplots(1, 1, figsize=(18, 6))
        ax = axes
        
        param_i = param_mcmc[i] 
        sample_t = np.array([samples_mcmc.transpose()[i]])
        sample_t = sample_t.transpose()
        if param_i=="dec_image" or param_i=="ra_image":
            param_i+="_"+str(i)
            if param_i=="dec_image":
                n_dec+=1
            else:
                n_ra+=1
        num_avrg = 100
        num_samples=len(sample_t[:,0])
        if int(num_samples/num_avrg)>1000:
            num_avrg= num_samples/960
        if j>= len(colors):
            j=0
        col=colors[j]
        my_plot_mcmc_behaviour(ax=axes,samples_mcmc=sample_t,param_mcmc=[param_i], num_average=num_avrg,col=col)
        del_ax = ax.plot([],[])
        del_ax[0].set_color("r")
        title_mcmc_behaviour = 'MCMC behaviour for '+param_i

        ax.set_title(title_mcmc_behaviour)
        plt.savefig(savefig_path+'MCMC_behaviour_'+param_i+'.png')
        plt.close()
        
        #MOD_BOUNDARIES
        if "ra_image" in param_i:
            min_bound,max_bound = find_bound(param_i, setting,num=n_ra-1)
        elif "dec_image" in param_i:
            min_bound,max_bound = find_bound(param_i, setting,num=n_dec-1)
        else:
            min_bound,max_bound = find_bound(param_i, setting)
            
            
        f, axes = plt.subplots(1, 1, figsize=(18, 6))
        ax = axes
        my_plot_mcmc_simil_behaviour(ax=axes,samples_mcmc=sample_t,param_mcmc=[param_i], 
                                     boundary=[min_bound,max_bound],
                                     num_average=num_avrg,col=col)
        del_ax = ax.plot([],[])
        del_ax[0].set_color("r")
        title_mcmc_behaviour = 'MCMC simil behaviour for '+param_i

        ax.set_title(title_mcmc_behaviour)
        plt.savefig(savefig_path+'simil_mcmc_bhv/MCMC_behaviour_'+param_i+'_simil.png')
        plt.close()
        j+=1
        


    # In[67]:


    # Time delay behaviour
    ########################
    from Utils.Dt_from_Df_reworked import Dt_XY
    param_mcmc_dt = ["Dt AB","Dt AC","Dt AD"]
    #D_fermat ordered as D_AB, D_AC, D_AD
    mcmc_Df = (np.array(mcmc_fermat.T[1:])-mcmc_fermat.T[0]).T
    H0_planck = 67.4 #km/s/Mpc
    mcmc_Dt = Dt_XY(Df_XY=mcmc_Df,H0=H0_planck,z_l=setting.z_lens,z_s=setting.z_source)
    samples_mcmc = mcmc_Dt
    param_mcmc   = param_mcmc_dt
    j=0
    for i in range(len(param_mcmc_dt)):  
        f, axes = plt.subplots(1, 1, figsize=(18, 6))
        ax = axes
        
        param_i = param_mcmc_dt[i] 
        sample_t = np.array([mcmc_Dt.transpose()[i]])
        sample_t = sample_t.transpose() 
        num_avrg = 100
        num_samples=len(sample_t[:,0])
        if int(num_samples/num_avrg)>1000:
            num_avrg= num_samples/960
        if j>= len(colors):
            j=0
        col=colors[j]
        my_plot_mcmc_behaviour(ax=axes,samples_mcmc=sample_t,param_mcmc=[param_i], num_average=num_avrg,col=col)
        del_ax = ax.plot([],[])
        del_ax[0].set_color("r")
        title_mcmc_behaviour = 'MCMC behaviour for '+param_i

        ax.set_title(title_mcmc_behaviour)
        plt.savefig(savefig_path+'MCMC_behaviour_'+param_i+'.png')
        plt.close()
        
        f, axes = plt.subplots(1, 1, figsize=(18, 6))
        ax = axes
        my_plot_mcmc_simil_behaviour(ax=axes,samples_mcmc=sample_t,param_mcmc=[param_i], num_average=num_avrg,col=col)
        del_ax = ax.plot([],[])
        del_ax[0].set_color("r")
        title_mcmc_behaviour = 'MCMC simil behaviour for '+param_i

        ax.set_title(title_mcmc_behaviour)
        plt.savefig(savefig_path+'simil_mcmc_bhv/MCMC_behaviour_'+param_i+'_simil.png')
        plt.close()


    # In[ ]:


    #MOD_LOG
    ########
    if logL_go:

            def my_plot_mcmc_logL(ax, mcmc_logL, num_average=100):
                """
                plots the MCMC behaviour and looks for convergence of the chain
                :param samples_mcmc: parameters sampled 2d numpy array
                :param param_mcmc: list of parameters
                :param dist_mcmc: log likelihood of the chain
                :param num_average: number of samples to average (should coincide with the number of samples in the emcee process)
                :return:
                """
                num_samples = len(mcmc_logL)
                if all_logL:
                    num_average = 1
                else:
                    num_average = int(num_average)
                n_points = int((num_samples - num_samples % num_average) / num_average)
                mcmc_logL_averaged = -np.max(mcmc_logL[:int(n_points * num_average)].reshape(n_points, num_average), axis=1)
                mcmc_logL_normed = (mcmc_logL_averaged - np.max(mcmc_logL_averaged)) / (np.max(mcmc_logL_averaged) - np.min(mcmc_logL_averaged))
                if all_logL:
                    ax.scatter(np.arange(0,len(mcmc_logL_normed)),mcmc_logL_normed,color='k',marker=".", alpha=0.001)
                else:
                    ax.plot(mcmc_logL_normed, color='k', linewidth=2)
                ax.set_title("Log Likelihood")
                ax.set_xlabel("MCMC Steps/"+str(num_average))
                return ax
            if np.inf in np.abs(mcmc_logL):
                print(sys.argv[0]+" - my_plot_mcmc_logL - WARNING: I am deleting part of the numbers for the LogL plot")
                print("n* of inf points: ",len(mcmc_logL[np.isinf(mcmc_logL)]) )
                inf_percentage = 100*np.round(len(mcmc_logL[np.isinf(mcmc_logL)])/len(mcmc_logL),4)
                print("percentage of inf points: ",inf_percentage,"%")
                if inf_percentage>50:
                    raise RuntimeError("inf point>50% of total points, something is wrong")
                from copy import deepcopy
                tmp_logL = deepcopy(mcmc_logL)
                tmp_logL [tmp_logL  == -np.inf] = max(tmp_logL)
                mcmc_logL[mcmc_logL == -np.inf] = min(tmp_logL)
                del tmp_logL
            f, axes = plt.subplots(1, 1, figsize=(18, 6))
            ax = axes
            num_avrg = 100
            num_samples=len(mcmc_logL)
            if int(num_samples/num_avrg)>1000:
                num_avrg= num_samples/960
            my_plot_mcmc_logL(axes,mcmc_logL,num_avrg)
            del_ax = ax.plot([],[])
            del_ax[0].set_color("r")
            plt.savefig(savefig_path+'/MCMC_bhv_LogL.png')
            plt.close()
            
            if all_logL:
                from corner import quantile
                #plot histogram comparing posterior distribution

                plt.figure(figsize=(12,9))    
                top = len(mcmc_logL)
                max_data = max(mcmc_logL)
                min_data = min(mcmc_logL)
                small_diff_min_max= (max_data-min_data)*.01
                n_bins= 3*int(round( 1 + 3.322*np.log10(top)))
                n_bins = np.linspace(min_data-small_diff_min_max,max_data+small_diff_min_max,n_bins+2 )

                mean = np.mean(mcmc_logL)
                sigma = np.std(mcmc_logL)
                val_min1sig, val, val_max1sig = quantile(mcmc_logL,q=[0.16, 0.5, 0.84])
                val_min3sig, val, val_max3sig = quantile(mcmc_logL,q=[0.0027,0.5,0.9973])

                sig_min1sig = np.abs(val_min1sig-val)
                sig_max1sig = val_max1sig - val

                sig_min3sig = np.abs(val_min3sig-val)
                sig_max3sig = val_max3sig - val


                count, _ , _ =plt.hist(mcmc_logL, bins=n_bins, density=True,alpha=0.2,color='k',
                         label=str(round(mean,2))+r"$^{+"+str(round(sig_max1sig,3))+"}_{-"+str(round(sig_min1sig,3))+"}"+\
                        "$ $^{+"+str(round(sig_max3sig,3))+"}_{-"+str(round(sig_min3sig,3))+"}$")
                y_max = max(count)*1.1
                plt.ylim(0,y_max)


                plt.vlines(val,0,y_max, color="b",label=r"Mean ",linestyles="dashed")
                plt.vlines(val_min1sig,0,y_max, color="r",label=r"1$\sigma$ region",linestyles="dashed")
                plt.vlines(val_max1sig,0,y_max, color="r",linestyles="dashed")
                plt.vlines(val_min3sig,0,y_max, color="k",label=r"3$\sigma$ region",linestyles="dashed")
                plt.vlines(val_max3sig,0,y_max, color="k",linestyles="dashed")

                plt.title("Distribution of logL")
                plt.xlabel("LogL []")
                plt.legend(loc="upper right")
                
                plt.savefig(savefig_path.replace("MCMC_bhv","")+"all_PD/PD_logL.png")
                plt.close() 


    # In[42]:


    success(sys.argv[0])

