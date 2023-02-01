#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import corner
import argparse
import numpy as np
import matplotlib.pyplot as plt

from tools import *

from Prior import mag_prior,mag_prior_ABC
from mag_remastered import labels
from pycs_get_res import get_combined_mag
from my_pycs_scripts.from_Dmag_to_FR import get_FR,get_sig_FR

def from_list_str_to_list(list_str_i):
    if list_str_i not in ["[","]",","]:
         return(float(list_str_i) )
         
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the prior of the mag ratio given the setting file")
    parser.add_argument("-PP", "--plot_prior", action="store_true", dest="plot_prior", default=False,
                        help="Plot the prior distribution for the mag ratio at the prior image positions")
    parser.add_argument("-ABCD", action="store_true", dest="ABCD", default=False,
                        help="Consider images A,B, C and also D (usually ignored due to low S/N)")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=5000,
                        help="Number of points for the MCMC prior")
    parser.add_argument("-magB", "--mag_boundaries",dest="mag_boundaries", default=False,action="store_true",
                            help="Boundaries of the mag ratio space")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args           = parser.parse_args()
    npoints        = args.npoints
    plot_prior     = args.plot_prior
    ABCD           = args.ABCD
    mag_boundaries = args.mag_boundaries
    if mag_boundaries is not False:
        # tmp
        mag_boundaries = [[-1 ,4],[ -3 , 1.5],[ -3, 3]]   # np.reshape([from_list_str_to_list(i)  for i in mag_boundaries],(3,2))
    else:
        mag_boundaries = None 
    settings       = args.SETTING_FILES
    ####################
    present_program(sys.argv[0])
    ####################
    backup_path   = "./backup_results/"
     
    prior_name = "mcmc_mag_prior"
        
    if not ABCD:
        labels[-1] = "$\mu_C$/$\mu_B$"
        prior_name = prior_name+"_ABC"
        lbl = ["AB","AC","BC"]
    else:
        lbl = ["AB","AC"]
        
    for sets in settings:
        print_setting(sets)
        sets = get_setting_module(sets,1)
        pycs_magcomb = get_combined_mag(sets.pycs_lensname+"_"+sets.pycs_dataname,main_dir_path=sets.pycs_path)
        FR = [get_FR(res) for res in pycs_magcomb.results]
        Rmag_ABC = []
        for i,lbl_i in enumerate(pycs_magcomb.labels):
            if lbl_i in lbl:
                Rmag_ABC.append(FR[i])
        if not ABCD:
            mag_prior_i = mag_prior_ABC(sets,npoints,save_mcmc=True,mag_boundaries=mag_boundaries,backup_path=backup_path,output_name=prior_name)
        else:    
            mag_prior_i = mag_prior(sets,npoints,save_mcmc=True,mag_boundaries=mag_boundaries,backup_path=backup_path,output_name=prior_name)
        if plot_prior:
            savefig_path = get_savefigpath(sets,backup_path)
            print("Expected FR:" , Rmag_ABC)
            corner.corner(np.abs(mag_prior_i),labels=labels,show_titles=True,truths=Rmag_ABC)
            plt.savefig(savefig_path+"/"+prior_name+".png")
            print("Saving plot ",savefig_path+"/"+prior_name+".png")
            plt.close()
    success(sys.argv[0])
