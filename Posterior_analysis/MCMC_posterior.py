#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# MOD_KS: I implement an analysis of the last 1/4 and second-to-last 
# 1/4 of the chain and do a KS statistic over it to see if it is compatible


# In[ ]:


# See old/MCMC_posterior.ipynb


# In[ ]:


import sys
import numpy as np
import json
from scipy import stats
import os
from datetime import datetime
import time 
import matplotlib.pyplot as plt
import pickle
from corner import quantile
import argparse
import importlib
from Utils.tools import *
from Utils.get_res import get_mcmc

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot the posterior and save the 3-sigma in a RDB file")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-e", "--evol",  dest="evol", default=False,action="store_true",
                        help="Plot the evolution of the posterior distribution")
    parser.add_argument("-n", "--norm",dest="norm", default=False, action="store_true",
                        help="Set the normalisation of the evolution distributions")
    parser.add_argument("-fc", "--full_corner", dest="full_corner", default=False,action="store_true",
                        help="Produces the standard corner plots")
    #MOD_KS
    parser.add_argument("-ks", "--KS", dest="ks", default=False,action="store_true",
                        help="Apply the KS statistic between the two last 1/4 of the posterior distribution")
    parser.add_argument('SETTING_FILE', help="setting file to consider")
    args = parser.parse_args()
    setting_name = get_setting_name(args.SETTING_FILE)
    setting      = get_setting_module(setting_name)
    backup_path  = "backup_results/"

    cut_mcmc = args.cut_mcmc
    full_corner = args.full_corner
    evol = args.evol
    norm = args.norm
    ks   = args.ks 
    
    ####################
    present_program(sys.argv[0])
    ####################
    kw_mcmc     = get_mcmc(setting_name,backup_path)

    samples     = kw_mcmc["mcmc_smpl"]
    param_names = kw_mcmc["mcmc_prm"]

    savefig_path = get_savefigpath(setting)
    #For the moment we skip the image position, 
    # see old/MCMC_posterior.ipynb for previous version

    # I also consider the standard corner plot
    if full_corner:
        import corner
        plot = corner.corner(samples, labels=param_names, show_titles=True)
        plot.savefig(savefig_path+"MCMC_posterior.png")
 
    str_param = setting.str_param

    all_param = []
    for i in param_names:
            if "ra_" in i or "dec_" in i:
                break
            else:
                all_param.append(i)
    col = ["r","b","g","y","k"]

    print_res = open(savefig_path+"res_3sigma.rdb","w")
    print_res.write("#3 Sigma results for "+setting_name+" \n")
    print_res.write("#################################\n\n")
    print_res.write("Parameter   Value   1-sigma(-)   1-sigma(+)   3-sigma(-)   3-sigma(+)\n")
    print_res.write("=========   =====   ==========   ==========   ==========   ==========\n")

    mkdir(savefig_path+"/all_PD")
    
    


    for i in range(len(all_param)):
        param = all_param[i]
        param_title, UOM = str_param(param)

        if param=="psi_ext":
            sample_i = np.transpose(samples)[i]*180/np.pi
        else:
            sample_i = np.transpose(samples)[i]  
        #plot histogram comparing posterior distribution
        plt.figure(figsize=(12,9))    
        top = len(sample_i)
        max_data = max(sample_i)
        min_data = min(sample_i)
        small_diff_min_max= (max_data-min_data)*.01
        n_bins= 3*int(round( 1 + 3.322*np.log10(top)))
        n_bins = np.linspace(min_data-small_diff_min_max,max_data+small_diff_min_max,n_bins+2 )

        mean = np.mean(sample_i)
        sigma = np.std(sample_i)
        val_min1sig, val, val_max1sig = quantile(sample_i,q=[0.16, 0.5, 0.84])
        val_min3sig, val, val_max3sig = quantile(sample_i,q=[0.0027,0.5,0.9973])

        sig_min1sig = np.abs(val_min1sig-val)
        sig_max1sig = val_max1sig - val

        sig_min3sig = np.abs(val_min3sig-val)
        sig_max3sig = val_max3sig - val


        count, _ , _ =plt.hist(sample_i, bins=n_bins, density=True,alpha=0.2,color=col[i%len(col)],
                 label=str(round(mean,2))+r"$^{+"+str(round(sig_max1sig,3))+"}_{-"+str(round(sig_min1sig,3))+"}"+\
                "$ $^{+"+str(round(sig_max3sig,3))+"}_{-"+str(round(sig_min3sig,3))+"}$")
        y_max = max(count)*1.1
        plt.ylim(0,y_max)

        #mean, 1-sigma and 3-sigma levels

        plt.vlines(val,0,y_max, color="b",label=r"Mean ",linestyles="dashed")
        plt.vlines(val_min1sig,0,y_max, color="r",label=r"1$\sigma$ region",linestyles="dashed")
        plt.vlines(val_max1sig,0,y_max, color="r",linestyles="dashed")
        plt.vlines(val_min3sig,0,y_max, color="k",label=r"3$\sigma$ region",linestyles="dashed")
        plt.vlines(val_max3sig,0,y_max, color="k",linestyles="dashed")

        plt.title("Posterior distribution of "+param_title)
        plt.xlabel(param_title+" "+UOM)
        plt.legend(loc="upper right")
        plt.savefig(savefig_path+"all_PD/PD_"+param+".png")
        #plt.show()
        plt.close()
        print_res.write(param+"    "+str(round(mean,2))+"    "+
                        str(round(sig_max1sig,3))+"    "+str(round(sig_min1sig,3))+
                "      "+str(round(sig_max3sig,3))+"    "+str(round(sig_min3sig,3))+"\n" )
    print_res.close()
    


    #####
    # evolution of PDF with steps
    #####
    if evol :
        import imageio
        #MOD_KS
        if ks:
            from scipy.stats import ks_2samp as ks
            print_ks = open(savefig_path+"KS_lstQrt.rdb","w")
            print_ks.write("# K-S statistic for each parameter between the two lasts 1/4 of the chain\n")
            print_ks.write("##########################################################################\n")
            print_ks.write("Param    KS Statistic      P-Value\n")
            print_ks.write("##################################\n")

        for i in range(len(all_param)):
            param = all_param[i]
            param_title, UOM = str_param(param)

            if param=="psi_ext":
                sample_i = samples[i]*180/np.pi
            else:
                sample_i = samples[i]  

            #Def boundaries:
            top = len(sample_i)
            max_data = max(sample_i)
            min_data = min(sample_i)
            small_diff_min_max= (max_data-min_data)*.01
            n_bins= 3*int(round( 1 + 3.322*np.log10(top)))
            n_bins = np.linspace(min_data-small_diff_min_max,max_data+small_diff_min_max,
                                 n_bins+2 )
            plt.figure(figsize=(18,11)) 
            count, _, _ =plt.hist(sample_i, bins=n_bins, density=norm)
            plt.close()
            y_rng = [0, max(count)*1.5]
            x_rng = [min(n_bins)*0.8,max(n_bins)*1.2]

            #split the sample in ~5 (n_frames)
            n_frames = 5

            ln_smpl = int(len(sample_i)/n_frames -len(sample_i)%n_frames)
            """
            i_smpl = 0
            sample_split= []
            while i_smpl<len(sample_i):
                j_smpl = i_smpl+ln_smpl
                if j_smpl+ln_smpl>len(sample_i):
                    j_smpl=len(sample_i)
                if pure_evol: 
                    sample_split.append(sample_i[0:j_smpl])
                else:
                    sample_split.append(sample_i[i_smpl:j_smpl])
                i_smpl=j_smpl
            """
            #sample_split = np.array([sample_i[j:j + n_frames] for j in range(0, len(sample_i), n_frames)]).T ->wrong
            sample_split=np.array([sample_i[j:j + ln_smpl] for j in range(0, len(sample_i), ln_smpl)])
            im_smpl = []
            for j in range(len(sample_split)):
                sample_j = sample_split[j]
                #plot histogram comparing posterior distribution
                fig = plt.figure(figsize=(18,11))    

                mean = np.mean(sample_j)
                sigma = np.std(sample_j)
                val_min1sig, val, val_max1sig = quantile(sample_j,q=[0.16, 0.5, 0.84])
                val_min3sig, val, val_max3sig = quantile(sample_j,q=[0.0027,0.5,0.9973])

                sig_min1sig = np.abs(val_min1sig-val)
                sig_max1sig = val_max1sig - val

                sig_min3sig = np.abs(val_min3sig-val)
                sig_max3sig = val_max3sig - val


                plt.hist(sample_j, bins=n_bins, density=norm,alpha=0.2,color=col[i%len(col)],
                         label=str(round(mean,3))+r"$^{+"+str(round(sig_max1sig,3))+"}_{-"+str(round(sig_min1sig,3))+"}"+\
                        "$ $^{+"+str(round(sig_max3sig,3))+"}_{-"+str(round(sig_min3sig,3))+"}$")

                plt.ylim(*y_rng)

                #mean, 1-sigma and 3-sigma levels

                plt.vlines(val,*y_rng, color="b",label=r"Mean ",linestyles="dashed")
                plt.vlines(val_min1sig,*y_rng, color="r",label=r"1$\sigma$ region",linestyles="dashed")
                plt.vlines(val_max1sig,*y_rng, color="r",linestyles="dashed")
                plt.vlines(val_min3sig,*y_rng, color="k",label=r"3$\sigma$ region",linestyles="dashed")
                plt.vlines(val_max3sig,*y_rng, color="k",linestyles="dashed")
                if j==0:
                    number="1st"
                elif j==1:
                    number="2nd"
                else:
                    number=str(j+1)+"th"
                plt.title("Posterior distribution of "+param_title+": "+number+" part")
                plt.xlabel(param_title+" "+UOM)
                plt.legend(loc="upper right")
                fig.canvas.draw()
                #convert the image to numpy array
                data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
                data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
                im_smpl.append(data)
                plt.close()
            # to compare the middle and the final shape, we add
            # 2 final frames(identicals, just to have more time to look at them)
            fig = plt.figure(figsize=(18,11))    
            plt.hist(sample_i, bins=n_bins, density=norm,alpha=0.15,color=col[i%len(col)],
                         label="Final distribution")
                     #str(round(mean,3))+r"$^{+"+str(round(sig_max1sig,3))+"}_{-"+str(round(sig_min1sig,3))+"}"+\
                     #   "$ $^{+"+str(round(sig_max3sig,3))+"}_{-"+str(round(sig_min3sig,3))+"}$")
            plt.vlines(np.mean(sample_i),*y_rng, color="b",label=r"Final mean",linestyles="dashed")

            half_run_sample =sample_i[0:int(len(sample_i)/2 - len(sample_i)%2)]
            half_run_mean = np.mean(half_run_sample)
            plt.hist(half_run_sample, bins=n_bins, density=norm,alpha=0.15,color="k",
                         label="Half-run distribution")
            plt.vlines(half_run_mean,*y_rng, color="k",label=r"Half-run mean",linestyles="dotted")

            plt.ylim(*y_rng)
            plt.title("Final and half-run Posterior distribution of  "+param_title)
            plt.xlabel(param_title+" "+UOM)
            plt.legend(loc="upper right")
            fig.canvas.draw()
            #convert the image to numpy array
            #data = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
            data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
            data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            im_smpl.append(data)
            im_smpl.append(data)
            plt.close()
            imageio.mimsave(savefig_path+"all_PD/evol_"+param+".gif",im_smpl,fps=1)

            #MOD_KS
            if ks:
                n_ks = 4
                ln_smpl_ks = int(len(sample_i)/n_ks -len(sample_i)%n_ks)
                sample_split_ks=np.array([sample_i[j:j + ln_smpl_ks] for j in range(0, len(sample_i), ln_smpl_ks)])
                ks_i = ks(sample_split_ks[-2],sample_split_ks[-1])
                print_ks.write(param+"  "+str(ks_i[0])+"  "+str(ks_i[1])+"\n")
    if ks:
        print_ks.close()
        
        
    success(sys.argv[0])


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


"""
# To implement better
import corner

for i in range(len(all_param)):
    for j in range(len(all_param)):
        if j!=i:            
            f, axes = plt.subplots(1, 1, figsize=(6, 6), sharex=False, sharey=False)
            ax = axes
            #ax.set_title("")
            x = sample[i]
            y = sample[j]
            param_x = all_param[i]
            param_title_x,UOM_x = str_param(param_x)
            param_y = all_param[j]
            param_title_y,UOM_y = str_param(param_y)
            
            ax.set_xlabel(param_title_x+" "+UOM_x)
            ax.set_ylabel(param_title_y+" "+UOM_y)
            corner.hist2d(x,y,ax=ax)
            #plt.savefig(savefig_path+"/corner_"+param_x+"_"+param_y) -> way too many (but it works)
            plt.show()
"""


# In[ ]:




