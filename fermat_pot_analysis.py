#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#  I want to correct the images position, which now is fixed, to be set to the one of the corresponding
#  mcmc step


# In[30]:


import sys
import numpy as np
import os
from lenstronomy.LensModel.lens_model import LensModel
import importlib
import corner
import json
import argparse

from get_res import *
from tools import *
from order_images import get_new_image_order

#labels_Fermat = ["$\phi_{Fermat}A$", "$\phi_{Fermat}B$","$\phi_{Fermat}C$" ,"$\phi_{Fermat}D$"]
#labels_Df     = ["$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} AD$"]
labels_Fermat = ["$\phi_A$", "$\phi_B$","$\phi_C$" ,"$\phi_D$"]
labels_Df     = ["$\Delta\phi_{AB}$", "$\Delta\phi_{AC}$","$\Delta\phi_{AD}$"]


def gen_mcmc_fermat(mcmc,param_mcmc,setting):
    # mcmc_prior shape = param, n_points
    # param_mcmc = list of name of params  
    # setting   = setting module
    
    if not check_if_CP(setting):
        lens_model_list = ['SIE']
    else:     
        print("WARNING: Considering the PEMD profile for the main lens")
        lens_model_list = ['PEMD']
    lens_model_list =[*lens_model_list,'SIS','SHEAR_GAMMA_PSI']
    
    lensModel = LensModel(lens_model_list=lens_model_list, 
                          z_lens=setting.z_lens,
                          z_source=setting.z_source)

    mcmc_fermat = []
    for i in range(len(mcmc)):
        kwargs_result_i = setting.produce_kwargs_result(mcmc,param_mcmc,i)
        x_image= kwargs_result_i["kwargs_ps"][0]["ra_image"][0]
        y_image= kwargs_result_i["kwargs_ps"][0]["dec_image"][0]
        fermat_potential = lensModel.fermat_potential(x_image,y_image,kwargs_result_i["kwargs_lens"])
        mcmc_fermat.append(fermat_potential.tolist())
    #I want to obtain the correct image order
    ########################################
    new_order = get_new_image_order(setting,mcmc)
    tmp_mcmc = np.array(mcmc_fermat).transpose()
    mcmc_fermat = np.transpose([tmp_mcmc[i] for i in new_order])
    return mcmc_fermat.tolist() # shape: (steps, f(i) )

def gen_mcmc_Df(mcmc,param_mcmc,setting):
    mcmc_fermat = gen_mcmc_fermat(mcmc,param_mcmc,setting)
    mcmc_Df = np.transpose(mcmc_fermat)[1:]-np.transpose(mcmc_fermat)[0] 
    mcmc_Df = mcmc_Df.T.tolist()  #shape: steps, D_AB, D_AC, D_AD, meaning Df_i - Df_A
    return mcmc_Df

"""
from corner import quantile

##########
# For DT #
##########
outputfile = [savefig_path+'read_results.data',
         savefig_path+'read_results_updated.data']

param_mcmc_dt = ["Dt_AB","Dt_AC","Dt_AD"]

import scipy.integrate as integrate
from astropy import units as u
from astropy import constants as const

Omega_m = 0.307       # omega matter
Omega_L = 1-Omega_m   # omega lambda 
Omega_k = 0.0         # omega curvature ->ignored from now on
H0_planck = 67.4*u.km/u.second/u.megaparsec # km/s/Mpc

def E(z,Omega_m0=Omega_m,Omega_L0=Omega_L):
    return (Omega_m0*((1+z)**3) + Omega_L0)**0.5  #non-dimensional

def k(z_l,z_s): #also = (1+z_d) * D_d D_d /D_ds
    zs    = integrate.quad(lambda x: (1./E(x)) ,0,z_s)[0]
    zl    = integrate.quad(lambda x: (1./E(x)) ,0,z_l)[0]
    zl_zs = integrate.quad(lambda x: (1./E(x)) ,z_l,z_s)[0]
    return zs*zl/zl_zs #non-dimensional

def delta_t(D_phi,z_l,z_s,H_0):
    D_phi *=(u.arcsecond.to("radian"))**2 # from arcsec^2 to rad^2
    H_0 = H_0.to_value("Hz")
    dt  = D_phi*k(z_l,z_s)/H_0 #sec
    dt *= u.second.to("day") #day
    return (dt)


print_res = open(savefig_path+"Single_fermat_potential.txt","w")
kwargs_fermat={}
kwargs_Df={}
print_res.write("We are considering the settings_file:" +str(setting_name))
print_res.write("\n#################################\n\n")

for i in range(len(labels_new)):
    val_min, val, val_max = quantile(mcmc_fermat[:,i],q=[0.16, 0.5, 0.84])
    sig_min = np.abs(val_min-val)
    sig_max = val_max - val
    sig_ref = np.min([sig_max,sig_min])
    n_exp = np.log10(1/sig_ref)
    fact=pow(10,int(n_exp)+2)
    if sig_min==sig_max:
        print_res.write(labels_new[i]+"  " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " +- "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    else:
        print_res.write(labels_new[i]+" " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " - "+str(np.trunc(np.array(sig_min)*fact)/fact)+\
                        " + "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")

    kwargs_fermat[labels_new[i]]=(val,sig_min,sig_max)
print_res.write("\n#################################\n")
print_res.write("\nFermat potential differences \n")
for i in range(len(labels_Df)):            
    val_min, val, val_max = quantile(np.array(mcmc_Df)[:,i],q=[0.16, 0.5, 0.84])
    sig_min = np.abs(val_min-val)
    sig_max = val_max - val
    sig_ref = np.min([sig_max,sig_min])
    n_exp = np.log10(1/sig_ref)
    fact=pow(10,int(n_exp)+2)
    if sig_min==sig_max:
        print_res.write(labels_Df[i]+"  " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " +- "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
    else:
        print_res.write(labels_Df[i]+" " +str(np.trunc(np.array(val)*fact)/fact)+\
                        " - "+str(np.trunc(np.array(sig_min)*fact)/fact)+\
                        " + "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")

    kwargs_Df[labels_Df[i]]=(val,sig_min,sig_max)
print_res.write("\n#################################\n")
print_res.close()
"""


# In[ ]:


if __name__=="__main__":
    #############################
    present_program(sys.argv[0])
    #############################
    parser = argparse.ArgumentParser(description="Produces the stacked MCMC results for the fermat potential \
                                     at the position of the images")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    
    args = parser.parse_args()

    setting_names =  args.SETTING_FILES.replace(".py","")
    cut_mcmc = int(args.cut_mcmc)
    backup_path="backup_results"

    for setting_name in setting_names:
        ############################################################
        #This should be the same for all settings
        if not check_if_CP(setting_name):
            lens_model_list = ['SIE']
        else:     
            print("WARNING: Considering the PEMD profile for the main lens")
            lens_model_list = ['PEMD']
        lens_model_list= [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']
        ############################################################

        savemcmc_path = get_savemcmcpath(setting_name)
        savefig_path  = get_savefigpath(setting_name)

        ############################################
        setting_module = get_setting_module(setting_name) 
        setting        = setting_module.setting()
        ############################################


        #MCMC sample
        samples_mcmc = get_mcmc_smpl(setting_name,backup_path)

        #parameters' name
        param_mcmc  = get_mcmc_prm(setting_name,backup_path)

        mcmc_fermat = gen_mcmc_fermat(samples_mcmc,param_mcmc,setting)
        mcmc_DfT    = np.transpose(mcmc_fermat)[1:]-np.transpose(mcmc_fermat)[0] 
        mcmc_Df     = mcmc_DfT.T
        
        #Save the mcmc in a file, NOTE: they are ordered A,B,C,D
        mcmc_file_name=savemcmc_path+setting_name.replace(".py","").replace("settings","mcmc_ordered_fermat")+".json"
        with open(mcmc_file_name, 'w+') as mcmc_file:
            json.dump(mcmc_fermat, mcmc_file)

        mcmc_fermat = np.array(mcmc_fermat[cut_mcmc:])
        mcmc_Df     = mcmc_Df[cut_mcmc:]

        plot = corner.corner(mcmc_fermat, labels=labels_Fermat, show_titles=True)
        plot.savefig(savefig_path+"Single_fermat_potential.png")


        plot = corner.corner(mcmc_Df, labels=labels_Df, show_titles=True)
        plot.savefig(savefig_path+"Single_Df.png")

    success(sys.argv[0])

