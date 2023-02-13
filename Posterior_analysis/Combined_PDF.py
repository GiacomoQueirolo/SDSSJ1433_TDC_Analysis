#!/usr/bin/env python
# coding: utf-8

# In[6]:


import numpy as np
import json
from scipy import stats
import os
from datetime import datetime
import time 
import matplotlib.pyplot as plt
import pickle
from corner import quantile
import sys
import pathlib as pth
import argparse
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
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()

    setting_name = [st.replace(".py","") for st in args.SETTING_FILES]
    cut_mcmc = int(args.cut_mcmc)
    inv_cut  = int(args.inv_cut)
    dir_name = args.dir_name
    backup_path = "backup_results"

    CP = []
    for st in setting_name:
        CP.append(check_if_CP(get_setting_module(st).setting()))
    CP = np.all(CP)

    if CP:
        print("WARNING: Considering the Change Profile -> PEMD")

    filters=[]
    for i in setting_name:
        filters.append(i.replace("settings_",""))

    if not CP:
        save_plot = "PDF_superposition_II"
    else:
        save_plot = "PDF_superposition_CP"
    save_plot = create_dir_name(setting_name,save_dir=save_plot,dir_name=dir_name,backup_path=backup_path)

    savemcmc_path=[]
    mcmc_file_name =[]
    param_file_name =[]
    for i in setting_name:
        path_i="./"+backup_path+"/"+i.replace("settings_","mcmc_")+"/"
        savemcmc_path.append(path_i)
        mcmc_file_name.append(path_i+i.replace("settings","mcmc_smpl")+".json")
        param_file_name.append(path_i+i.replace("settings","mcmc_prm")+".dat")


    for st in setting_name:
        setting_dir = find_setting_path(st)
        os.system("cp "+setting_dir+"/"+str(st)+".py "+str(save_plot)+"/.")


    samples = []
    for i in mcmc_file_name:
        with open(i, 'r') as mcmc_file_i:
            samples.append(np.array(json.load(mcmc_file_i)).T)


    param_names=[]
    for p in param_file_name:
        with open(p,"r") as param_file:
            pn=(param_file.readlines())
        for i in range(len(pn)):
            pn[i]=pn[i].replace(",\n","")
        param_names.append(pn)

    # for fermat potentials
    fermat_mcmc=[get_mcmc_fermat(st,backup_path) for st in setting_name]


    samples_Df = []
    for i in fermat_mcmc:
        with open(i, 'r') as mcmc_file_i:
            mcmc_i  =np.array(json.load(mcmc_file_i)).T
        mcmc_Dfi=np.array(mcmc_i[1:]-mcmc_i[0])
        samples_Df.append(mcmc_Dfi)
            
    for i in range(len(samples)):
        samples[i] = np.array(samples[i].tolist() + samples_Df[i].tolist())
        param_names[i]= param_names[i]+["$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} AD$"]

        
    # For mag ratio

    samples_mag = []
    for i in range(len(setting_name)):
        mag_mcmc_file = savemcmc_path[i]+setting_name[i].replace("settings","mcmc_mag_rt")+".json"

        if not os.path.isfile(mag_mcmc_file):
            raise FileNotFoundError("Cannot find "+mag_mcmc_file)
            
        with open(mag_mcmc_file,"r") as f:
            mcmc_mag_ratio = json.load(f)
        samples_mag.append(np.transpose(mcmc_mag_ratio))

        
    for i in range(len(samples)):
        samples[i] = np.array(samples[i].tolist() + samples_mag[i].tolist())
        param_names[i]= param_names[i]+["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_D$/$\mu_A$"]
        
    #For the moment we skip the image position -> see old/Combined_PDF.py #1


    # In[ ]:


    param_to_compare= ['theta_E_lens0', 'e1_lens0', 'e2_lens0','theta_E_lens1', 'gamma_ext_lens2','psi_ext_lens2',\
                       "$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} AD$",\
                       "$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_D$/$\mu_A$"]
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
        


    # In[ ]:


    if cut_mcmc:
        print("\nCUT: Ignoring the first ",cut_mcmc," steps of the MCMC sample\n")
        samples = [samples[i].T[cut_mcmc:].T for i in range(len(samples))]
    if inv_cut:
        print("\nINVERSE CUT: Considering only the last ",inv_cut," steps of the MCMC sample\n")
        samples = [samples[i].T[-inv_cut:].T for i in range(len(samples))]
        


    # In[2]:


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
                UOM=r'[\"$^2$]'
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
        #n_bins = math.ceil((max_data - min_data)/w)
        n_bins = np.linspace(min_data-small_diff_min_max,max_data+small_diff_min_max,n_bins+2)
        for s in range(len(sample_comp)):
            sig_min,mean,sig_max = quantile(sample_comp[s],q=[0.16,0.5,0.84])
            sig_min,sig_max = np.array([mean-sig_min]), np.array([sig_max-mean])
            sigma = np.mean([sig_min,sig_max])        
            count, bins, _ = plt.hist(sample_comp[s], bins=n_bins, density=True,alpha=0.2,color=col[s],label=filters[s]+":"+str(round(mean,2))+r"$\pm$"+str(round(sigma,3)) )
            plt.errorbar(mean,max(count)/2.,yerr=None,xerr=[sig_min,sig_max],fmt=str(col[s]+"+"))
            plt.scatter(mean,max(count)/2.,c=col[s],marker=".")

        plt.title("Posterior distribution of "+param_title+" for "+lf)
        plt.xlabel(param_title+" "+UOM)
        plt.legend(loc="upper right")
        if "mu" in param:
            param = "mu_"+param.replace("\\","").replace("/","").replace("$","").replace("mu_","")[::-1]
        plt.savefig(str(save_plot)+"/SPD_"+param+".png")
        #plt.show()
        plt.close() 

    # In[ ]:


    # see old/Combined_PDF.py, section 2: obsolete way to combine the PDF


    # In[ ]:


    print("Result directory:", str(save_plot))

    success(sys.argv[0])

