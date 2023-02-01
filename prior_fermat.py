#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import corner
import argparse
import matplotlib.pyplot as plt

from tools import *
from Prior import Df_prior,Df_prior_ABC


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the prior of the difference of Fermat potential differences given the setting file")
    parser.add_argument("-PP", "--plot_prior", action="store_true", dest="plot_prior", default=False,
                        help="Plot the prior distribution for the Fermat potential difference at the prior image positions")
    parser.add_argument("-ABCD", action="store_true", dest="ABCD", default=False,
                        help="Consider images A,B, C and also D (usually ignored due to low S/N)")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=1000,
                        help="Number of points for the MCMC prior")
    parser.add_argument("-DfB", "--Df_boundaries",dest="Df_boundaries", default=None,
                        help="Boundaries of the DFermat potential space")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    npoints       = args.npoints
    plot_prior    = args.plot_prior
    ABCD          = args.ABCD
    Df_boundaries = args.Df_boundaries
    settings      = args.SETTING_FILES
    ####################
    present_program(sys.argv[0])
    ####################
    backup_path   = "./backup_results/"
    
    labels = ["$\Delta\phi_{AB}$","$\Delta\phi_{AC}$","$\Delta\phi_{AD}$"]
    prior_name = "mcmc_Df_prior"
    if not ABCD:
        labels[-1] = "$\Delta\phi_{BC}$"
        prior_name = prior_name+"_ABC"
    for sets in settings:
        print_setting(sets)
        if not ABCD:
            Df_prior_i = Df_prior_ABC(sets.replace(".py",""),npoints,save_mcmc=True,Df_boundaries=Df_boundaries,backup_path=backup_path,output_name=prior_name)
        else:    
            Df_prior_i = Df_prior(sets.replace(".py",""),npoints,save_mcmc=True,Df_boundaries=Df_boundaries,backup_path=backup_path,output_name=prior_name)
        if plot_prior:
            savefig_path = get_savefigpath(sets,backup_path)
            corner.corner(Df_prior_i,labels=labels,show_titles=True)
            plt.savefig(savefig_path+"/"+prior_name+".png")
            plt.close()
    success(sys.argv[0])

