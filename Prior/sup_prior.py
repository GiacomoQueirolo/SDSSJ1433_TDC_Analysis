#!/usr/bin/env python
# coding: utf-8

import sys
import argparse

from Utils.tools import *
from Prior.Prior import Prior,load_Prior
from Compare_Res.super_corner import plot_sup_corner_Df,plt_sup_corner_lnsprm

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Produces the prior of the lens parameters and the difference of Fermat potential \
                                     given the setting files and overplot them to show that they are similar enough")
    parser.add_argument("-AD", dest="AD", default=False,action="store_true",
                        help="Consider AD couple instead of BC")
    parser.add_argument("-DfB", "--Df_boundaries",dest="Df_boundaries", default=None,
                        help="Boundaries of the DFermat potential space")
    parser.add_argument("-sl","--simple_legend", action="store_true", dest="simple_legend", default=False,
                        help="Use simplified legend in the plot (only filters names)")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=1000,
                        help="Number of points for the MCMC prior")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default="",
                    help="Directory name where to save the plot")
    parser.add_argument("-owP","--overwrite_Prior", action="store_true", dest="overwrite_Prior", default=False,
                        help="If same prior is present (same N_prior), recalculate it overwrite it ")

    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    npoints       = args.npoints
    AD            = args.AD
    Df_boundaries = args.Df_boundaries
    settings      = args.SETTING_FILES
    simple_legend = args.simple_legend
    overwrite_Prior = bool(args.overwrite_Prior)
    dir_name      = str(args.dir_name)
    if dir_name!="" and dir_name[0]!="_":
        dir_name="_"+dir_name
    ####################
    present_program(sys.argv[0])
    ####################
    backup_path = "./backup_results/"
    save_dir    = create_dir_name(settings,save_dir="PDF_Sup_Corner/",dir_name="/prior"+dir_name,backup_path=backup_path)   
    save_dir_Df = create_dir_name(settings,save_dir=save_dir,dir_name="/Df",backup_path=".")
    save_log_command(save_dir)
    if AD:
        from Posterior_analysis.fermat_pot_analysis import labels_Df_AD as labels
    else:
        from Posterior_analysis.fermat_pot_analysis import labels_Df_BC as labels

    smpl_priors = []
    Df_priors   = []
    for sets in settings:
        prior_name = f"{save_dir}/prior_obj_{strip_setting_name(sets,filter=True)}.dll"
        prior = Prior(sets,Nsample=npoints,BC=not AD)
        prior = load_Prior(prior_path=prior_name,prior=prior,verbose=False,compute_Df=True,overwrite=overwrite_Prior)
        smpl_priors.append(prior.get_sample())
        Df_priors.append(prior.get_Df_sample())
    plt_sup_corner_lnsprm(settings,smpl_mcmc=smpl_priors,prm_mcmc=None,stnd_lnsprm=True,savefig_dir=save_dir,cut_mcmc=None,simple_legend=simple_legend)
    plot_sup_corner_Df(settings,Df_mcmc=Df_priors,savefig_dir=save_dir_Df,name_pdf="Prior_Df",cut_mcmc=None,simple_legend=simple_legend,BC=not AD)
    success(sys.argv[0])
