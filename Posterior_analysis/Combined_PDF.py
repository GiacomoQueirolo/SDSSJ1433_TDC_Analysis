#!/usr/bin/env python
# coding: utf-8

# In[6]:


import os,sys 
import argparse
import numpy as np
from corner import quantile
import matplotlib.pyplot as plt

from Utils.tools import *
from Utils.get_res import *
from Utils.get_param_title import newProfile

if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################

    parser = argparse.ArgumentParser(description="Plot the posterior distribution of the results in the single filter directories and the superpositions PD\
              saved in backup_results/PDF_superpostion_II/. (different if CP true)")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-ic", "--inverse_cut", type=int, dest="inv_cut", default=0,
                        help="\"inverse cut\": consider only the last <IC> steps of the mcmc")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                        help="Directory name where to save the superposed posteriors")
    parser.add_argument("-BC", dest="BC", default=False,action="store_true",
                        help="Consider BC couple instead of AD")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    BC   = args.BC
    AD_or_BC     = "BC" if BC else "AD"
    muAD_or_muBC = "$\mu_C$/$\mu_B$" if BC else "$\mu_D$/$\mu_A$" 

    setting_names = [get_setting_name(st) for st in args.SETTING_FILES]
    if setting_names==[]:
        raise ValueError("Must at least have one setting file")
    settings = get_setting_module(setting_names,1)
    cut_mcmc = int(args.cut_mcmc)
    inv_cut  = int(args.inv_cut)
    dir_name = args.dir_name
    backup_path = "backup_results"

    CP = np.all([check_if_CP(st) for st in settings])

    if CP:
        print("WARNING: Considering the Change Profile -> PEMD")

    filters=[get_filter(sett_i) for sett_i in settings]
        
    if not CP:
        save_plot = "PDF_superposition_II"
    else:
        save_plot = "PDF_superposition_CP"
        
    save_plot = create_dir_name(setting_names,save_dir=save_plot,dir_name=dir_name,backup_path=backup_path)
    
    samples     = []
    param_names = []
    # fermat potentials
    prm_dfi     = ["$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} "+AD_or_BC+"$"]
    #  mag ratio
    prm_mri     = ["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$",muAD_or_muBC]
    for sett in settings:
        smpl_i   = get_mcmc_smpl(sett,backup_path).T.tolist()
        smpl_dfi = np.transpose(get_mcmc_Df(sett,backup_path,noD=BC)).tolist()
        smpl_mri = get_mcmc_mag(sett,backup_path).T.tolist()
        samples.append(smpl_i+smpl_dfi+smpl_mri)
        prm_i = get_mcmc_prm(sett)
        param_names.append(prm_i+prm_dfi+prm_mri)
        
    #For the moment we skip the image position -> see old/Combined_PDF.py #1
    
    param_to_compare= ['theta_E_lens0', 'e1_lens0', 'e2_lens0','theta_E_lens1', 'gamma_ext_lens2','psi_ext_lens2',*prm_dfi,*prm_mri]
        
    if CP:
        param_to_compare = [param_to_compare[0],"gamma_lens0",*param_to_compare[1:] ]

    pair=[]
    for par in param_to_compare:
        pair_i=[]
        for i in range(len(param_names)):
            for j in range(len(param_names[i])):
                if par == param_names[i][j]:
                    pair_i.append(j)
                    break
        pair.append(pair_i)
        
        
    if cut_mcmc:
        print("\nCUT: Ignoring the first ",cut_mcmc," steps of the MCMC sample\n")
        samples = [samples[i].T[cut_mcmc:].T for i in range(len(samples))]
    if inv_cut:
        print("\nINVERSE CUT: Considering only the last ",inv_cut," steps of the MCMC sample\n")
        samples = [samples[i].T[-inv_cut:].T for i in range(len(samples))]
        
    lf = filters[0]+", "
    for i in filters[1:-1]:
        lf+=i+ ", "
    lf+= "and "+ filters[-1]

    col = ["r","b","g","y","k"]
  

    for i in range(len(pair)):
        param = param_names[0][pair[i][0]]
        if not "Fermat" in param and not "\mu" in param:
            param_title,UOM = newProfile(param)
        else:
            param_title=param
            if "Fermat" in param: 
                UOM=r'[arcsec $^2$]'
            else: 
                UOM="[]"
        sample_comp = []
        for j in range(len(samples)):
            if param=="psi_ext":
                sample_comp.append(samples[j][pair[i][j]]*180/np.pi)
            else:
                sample_comp.append(samples[j][pair[i][j]])
        #plot histogram comparing posterior distribution
        plt.figure(figsize=(12,9))    
        top = min([len(k) for k in sample_comp])
        max_data = max([max(k) for k in sample_comp])
        min_data = min([min(k) for k in sample_comp])
        small_diff_min_max= (max_data-min_data)*.01
        n_bins= 3*int(round( 1 + 3.322*np.log10(top)))
        n_bins = np.linspace(min_data-small_diff_min_max,max_data+small_diff_min_max,n_bins+2)
        for s in range(len(sample_comp)):
            sig_min,mean,sig_max = quantile(sample_comp[s],q=[0.16,0.5,0.84])
            sig_min,sig_max = np.array([mean-sig_min]), np.array([sig_max-mean])
            sigma = np.mean([sig_min,sig_max])
            count, bins, _ = plt.hist(sample_comp[s], bins=n_bins, density=True,stacked=True,alpha=0.2,color=col[s],label=filters[s]+":"+str(round(mean,2))+r"$\pm$"+str(round(sigma,3)) )
            plt.errorbar(mean,max(count)/2.,yerr=None,xerr=[sig_min,sig_max],fmt=str(col[s]+"+"))
            plt.scatter(mean,max(count)/2.,c=col[s],marker=".")

        plt.title("Posterior distribution of "+param_title+" for "+lf)
        plt.xlabel(param_title+" "+UOM)
        plt.legend(loc="upper right")
        if "mu" in param:
            param = "mu_"+param.replace("\\","").replace("/","").replace("$","").replace("mu_","")[::-1]
        name_fig = str(save_plot)+"/SPD_"+param+".png"
        print("Saved plot in "+name_fig)
        plt.savefig(name_fig)
        plt.close() 

    print("Result directory:", str(save_plot))
    success(sys.argv[0])
