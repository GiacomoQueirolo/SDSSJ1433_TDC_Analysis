#!/usr/bin/env python
# coding: utf-8
import argparse

from Utils.tools import *
#from Prior.Prior import Df_prior,Df_prior_ABC
from Prior.Prior import Prior
from Compare_Res.super_corner import plot_sup_corner_Df
from Posterior_analysis.fermat_pot_analysis import labels_Df_AD,labels_Df_BC

raise NotImplementedError("Discontinued, use sup_prior.py instead")
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the prior of the difference of Fermat potential differences given the setting file")
    #parser.add_argument("-PP", "--plot_prior", action="store_true", dest="plot_prior", default=False,
    #                    help="Plot the prior distribution for the Fermat potential difference at the prior image positions")
    parser.add_argument("-AD", dest="AD", default=False,action="store_true",
                        help="Consider AD couple instead of BC")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=1000,
                        help="Number of points for the MCMC prior")
    #parser.add_argument("-DfB", "--Df_boundaries",dest="Df_boundaries", default=None,
    #                    help="Boundaries of the DFermat potential space")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                        help="Directory name where to save the prior plot")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    npoints       = args.npoints
    #plot_prior    = args.plot_prior
    BC            = not args.AD
    #Df_boundaries = args.Df_boundaries
    settings      = args.SETTING_FILES
    dir_name      = args.dir_name 
    ####################
    present_program(sys.argv[0])
    ####################
    backup_path   = "./backup_results/"
    savefig_dir   = create_dir_name(settings,save_dir="Prior_Df",dir_name=dir_name,backup_path=backup_path,copy_settings=True)
    save_log_command(savefig_dir)

    labels = labels_Df_AD 
    prior_name = "mcmc_Df_prior"
    if BC:
        labels = labels_Df_BC
        prior_name = prior_name+"_ABC"
    Dfs = []
    for sets in settings:
        print_setting(sets)
        prior_sett = Prior(setting=sets,BC=BC,Nsample=npoints)
        Df_prior_i = prior_sett.get_Df_sample()
        Dfs.append(Df_prior_i)
        #if BC:
        #    Df_prior_i = Df_prior_ABC(sets.replace(".py",""),npoints,save_mcmc=True,Df_boundaries=Df_boundaries,backup_path=backup_path,output_name=prior_name)
        #else:    
        #    Df_prior_i = Df_prior(sets.replace(".py",""),npoints,save_mcmc=True,Df_boundaries=Df_boundaries,backup_path=backup_path,output_name=prior_name)
        """if plot_prior:
            savefig_path = get_savefigpath(sets,backup_path)
            corner.corner(Df_prior_i,labels=labels,show_titles=True)
            plt.savefig(savefig_path+"/"+prior_name+".pdf")
            plt.close()"""
    
    plot_sup_corner_Df(setting_list=settings,Df_mcmc=Dfs,BC=BC,simple_legend=True,savefig_dir=savefig_dir,name_pdf=prior_name)
    success(sys.argv[0])
