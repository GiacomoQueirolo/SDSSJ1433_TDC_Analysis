#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import corner
import argparse
import matplotlib.pyplot as plt

from tools import *
from Prior import Df_prior,Df_prior_ABC
from super_corner import plot_sup_corner_Df

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the prior of the difference of Fermat potential differences given the setting files and overplot them to show that they are similar enough")

    parser.add_argument("-ABCD", action="store_true", dest="ABCD", default=False,
                        help="Consider images A,B, C and also D (usually ignored due to low S/N)")
    parser.add_argument("-DfB", "--Df_boundaries",dest="Df_boundaries", default=None,
                        help="Boundaries of the DFermat potential space")
    parser.add_argument("-sl","--simple_legend", action="store_true", dest="simple_legend", default=False,
                        help="Use simplified legend in the plot (only filters names)")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=1000,
                        help="Number of points for the MCMC prior")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default="",
                    help="Directory name where to save the plot")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    npoints       = args.npoints
    ABCD          = args.ABCD
    Df_boundaries = args.Df_boundaries
    settings      = args.SETTING_FILES
    simple_legend = args.simple_legend
    dir_name      = args.dir_name
    ####################
    present_program(sys.argv[0])
    ####################
    backup_path = "./backup_results/"
    save_dir    = create_dir_name(settings,save_dir="PDF_Sup_Corner/",dir_name="/prior_Df"+dir_name,backup_path=backup_path)    
    labels      = ["$\Delta\phi_{AB}$","$\Delta\phi_{AC}$","$\Delta\phi_{AD}$"]
    if not ABCD:
        labels[-1] = "$\Delta\phi_{BC}$"
    priors_Df = []
    for sets in settings:
        print_setting(sets)
        if not ABCD:
            Df_prior_i,_ = Df_prior_ABC(sets.replace(".py",""),npoints,save_mcmc=False,Df_boundaries=Df_boundaries)
        else:    
            Df_prior_i,_ = Df_prior(sets.replace(".py",""),npoints,save_mcmc=False,Df_boundaries=Df_boundaries)
        priors_Df.append(Df_prior_i)
    plot_sup_corner_Df(settings,fermat_mcmc=priors_Df,savefig_dir=save_dir,name_pdf="Prior_Df",cut_mcmc=None,simple_legend=simple_legend,param_names=labels,already_BC=True)
    success(sys.argv[0])

