#!/usr/bin/env python
# coding: utf-8

# Copy from Multiply_Df_PDF
# now considering fitting the Prior in the fermat pot space

import sys,os
import argparse
import numpy as np
import json,pickle,dill
from corner import corner
from lenstronomy.Util.sampling_util import sample_ball
from lenstronomy.Sampling.Pool.multiprocessing import MultiPool

#my libs
from Utils.tools import *
from Utils.get_res import *
from Prior.Prior import Prior
from Utils.statistical_tools import get_bins_volume
from Utils.combinedsetting_class import combined_setting
from Posterior_analysis.fermat_pot_analysis import get_mcmc_Df,labels_Df_AD,labels_Df_BC 

def create_combined_setting(settings,save_dir,cmb_sett_name=None,comment="",filters=None,BC=True):
    z_lens   = [s.z_lens   for s in settings ]
    z_source = [s.z_source for s in settings ]
    if filters is None:
        filters       = [get_filter(st) for st in settings]
    if all(zl==z_lens[0] for zl in z_lens) and all(zs==z_source[0] for zs in z_source):
        z_lens   = z_lens[0]
        z_source = z_source[0]
    else:
        raise RuntimeError("Not all settings have same redshift for lens and source. Something is off")
   
    CombSett      = combined_setting(comment,z_source,z_lens,filters,setting_names,savedir=save_dir,BC=BC)
    cmb_sett_name = CombSett.gen_cmb_sett_name(cmb_sett_name = cmb_sett_name)
    with open(f"combined_settings/{cmb_sett_name}.dll","wb") as f:
        dill.dump(CombSett,f)
    try:
        os.symlink(f"combined_settings/{cmb_sett_name}.dll",f"{save_dir}/{cmb_sett_name}.dll")
    except FileExistsError:
        if verbose:
            print(f"{cmb_sett_name} existed in {save_dir}, overwriting link")
        os.remove(f"{save_dir}/{cmb_sett_name}.dll")
        os.symlink(f"combined_settings/{cmb_sett_name}.dll",f"{save_dir}/{cmb_sett_name}.dll")
    return CombSett


if __name__=="__main__":


    parser = argparse.ArgumentParser(description="Plot the multiplied posterior distribution of the Fermat potential difference from the given filters",
                        formatter_class=CustomFormatter,
                        usage='pokus --help',)
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="Cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                        help="Directory name where to save the multiplied posteriors")
    parser.add_argument("-cmbn","--combined_setting_name",type=str,dest="cmb_sett_name", default=None,
                        help="Name for the combined setting file")
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=40,
                        help="Number of bins per dimension (be careful with it! too many bins can be catastrophic)")
    parser.add_argument("-ms","--mcmc_steps",type=int, dest="mcmc_steps", default=1000,
                        help="Number of steps for the MCMC sampling and plot")
    parser.add_argument("-mp","--mcmc_prior",type=int, dest="mcmc_prior", default=1000,
                        help="Number of steps for the MCMC sampling of the Priors")
    parser.add_argument("-KDE", action="store_true", dest="KDE", default=False,
                        help="Use KDE (Kernel Density Estimator) instead of histograms (WARNING:Very slow for high number of points and/or bins)")
    parser.add_argument("-mcmc","--MCMC", action="store_true", dest="mcmc", default=False,
                        help="Also do the MCMC integration of the posterior")  
    parser.add_argument("-fcNp","--factNprior",type=int, dest="factNprior", default=1,
                        help="Scaling factor for (N_prior), N_prior=Number of sampled prior points used for the fitting")
    parser.add_argument("-owP","--overwrite_Prior", action="store_true", dest="overwrite_Prior", default=False,
                        help="If same prior is present (same N_prior), recalculate it overwrite it ")
    parser.add_argument("-AD", action="store_true", dest="AD", default=False,
                        help="Consider AD couple instead of BC")
    parser.add_argument("-v","--verbose", action="store_true", dest="verbose", default=False,
                        help="Verbose")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

    args = parser.parse_args()
    cut_mcmc = int(args.cut_mcmc)
    dir_name = args.dir_name
    KDE = bool(args.KDE)
    nbins = int(args.nbins)
    mcmc  = bool(args.mcmc)
    mcmc_steps  = int(args.mcmc_steps)
    mcmc_prior  = int(args.mcmc_prior) 
    verbose     = bool(args.verbose)
    factNprior  = int(args.factNprior)
    setting_names =  args.SETTING_FILES  
    settings      = get_setting_module(setting_names,True)
    cmb_sett_name = str(args.cmb_sett_name)
    overwrite_Prior = bool(args.overwrite_Prior)
    BC              = not bool(args.AD)
    ########################################
    present_program(sys.argv[0])
    ########################################

    backup_path  = "backup_results"
    main_savedir = "PDF_multiplication_ABC"
    if not BC:
        main_savedir+="D"
    filters       = [get_filter(st) for st in settings]
    savemcmc_path = [get_savemcmcpath(st) for st in  settings]
    save_dir      = create_dir_name(setting_names,save_dir=main_savedir,dir_name=dir_name,backup_path=backup_path,copy_settings=True)

    save_log_command(save_dir)
    # for fermat potentials    
    samples = []
    param_names = labels_Df_AD
    if BC:
        param_names = labels_Df_BC
        print("WARNING: Given the low S/N of image D, I will discard here and instead consider Delta BC")
        
    for st in settings:
        Dfi = get_mcmc_Df(st,backup_path=backup_path,noD=BC)
        cut_mcmc_scaled = int(len(Dfi)*cut_mcmc/1000)
        DfiT = np.transpose(Dfi[cut_mcmc_scaled:])
        samples.append(DfiT) #shape: 3, len(mcmc)
    
    
    ##################################################################################
    # Prior #
    name_prior = "kw_Prior_"+"_".join(filters)+".pkl"
    path_prior = backup_path+"/"+main_savedir+"/"+name_prior
    Nsample_prior = int(factNprior*1e4)
    prior_name = f"{save_dir}/prior_obj.dll"
    prior = Prior(settings[0],Nsample=Nsample_prior,BC=BC)
    if not overwrite_Prior and os.path.isfile(prior_name):
        with open(prior_name,"rb") as f:
            loaded_prior = dill.load(f)
        if prior==loaded_prior:
            prior    = loaded_prior
            prior_df = loaded_prior.get_Df_sample()
        else:
            print("Prior found, but it's different from expected")
            if loaded_prior.lens_prior!=prior.lens_prior:
                raise RuntimeWarning(f"Previous prior is structurally different.")
            if loaded_prior.Nsample>prior.Nsample:
                raise RuntimeWarning(f"Previous prior {prior_name} had larger sample: {loaded_prior.Nsample} vs {Nsample_prior}. No default way to handle this. If so you have to delete it by hand")
            else:
                if verbose:
                    print(f"It has lower sample, moving it to old_prior_obj.dll")
            os.system(f"mv {prior_name} {save_dir}/old_prior_obj.dll")
            prior_df = prior.get_Df_sample()
            with open(prior_name,"wb") as f:
                dill.dump(prior,f)
    else:
        prior_df = prior.get_Df_sample()
        with open(prior_name,"wb") as f:
            dill.dump(prior,f)
    ##################################################################################


    if KDE:
        from Multiply_PDF import Multiply_PDF_KDE
        npoints = nbins  # To TEST
        #Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Priors=Priors_Df,savedir=save_dir)
        #Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Prior_LR=Prior_Df_LRes,Prior_HR=Prior_Df_HRes,savedir=save_dir)
        raise RuntimeWarning("Not implemented yet")
        Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Prior=Prior,savedir=save_dir)
    else:
        from Multiply_PDF import Multiply_PDF_HIST_fitPrior
        #Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Priors=Priors_Df,savedir=save_dir)
        #Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Prior_LR=Prior_Df_LRes,Prior_HR=Prior_Df_HRes,savedir=save_dir)
        Combined_PDF,Combined_bins = Multiply_PDF_HIST_fitPrior(samples,nbins,Prior=prior,savedir=save_dir,labels=param_names,verbose=verbose)

    # The Combined_PDF must be re-normalised
    if KDE:
        #Combined_PDF = Normalise_KDE(Combined_PDF,Positions_KDE)
        Combined_PDF /= np.sum(Combined_PDF)
        # still unclear how to do it in KDE
    else:
        Combined_PDF /= np.sum(Combined_PDF*get_bins_volume(Combined_bins))


    with open(str(save_dir)+"/Combined_PDF"+["_KDE" if KDE else ""][0]+".pkl","wb") as f:
        pickle.dump(Combined_PDF,f)
        
    if KDE:
        with open(str(save_dir)+"/Combined_PDF_KDE_positions.pkl","wb") as f:
            pickle.dump(Positions_KDE,f)
    else:
        with open(str(save_dir)+"/Combined_PDF_bins.pkl","wb") as f:
            pickle.dump(Combined_bins,f)

    #####################################################################################
    comment  = "Product of posteriors done by "+str(sys.argv[0])
    CombSett  = create_combined_setting(settings=settings,cmb_sett_name=cmb_sett_name,comment=comment,filters=filters,save_dir=save_dir,BC=BC)
    #####################################################################################

    # We need to sample it for the plot 
    if not KDE and mcmc:
        from Utils.statistical_tools import *
        from emcee import EnsembleSampler
        from multiprocessing import cpu_count
        mcmc_init_pos,mcmc_sigma = estimate_for_mcmc(Combined_PDF,Combined_bins)
        #mcmc_sampling = sampler(mcmc_init_pos,prob=Combined_PDF,bins=Combined_bins,mcmc_simga=mcmc_sigma,mcmc_steps=int(1e6))
        #mcmc_chain = mcmc_sampling[0]
        def logP(pos):
            prob_at_pos = get_prob_at_pos(pos,Combined_PDF,Combined_bins)
            if prob_at_pos==0:
                return -np.inf
            return np.log(prob_at_pos)

        nprocesses    = cpu_count()-1
        nwalkers      = 42
        init_sample   = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
        pool          = MultiPool(processes=nprocesses) 
        emcee_mcmc    = EnsembleSampler(nwalkers,len(np.shape(Combined_PDF)), logP,pool=pool)
        initial_state = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
        emcee_mcmc.run_mcmc(initial_state= initial_state,nsteps=mcmc_steps)
        mcmc_chain = emcee_mcmc.get_chain(discard=0, thin=1, flat=True)
        with open(str(save_dir)+"/mcmc_Df_chain.json","w") as f:
            json.dump(np.array(mcmc_chain).tolist(),f)
        plot = corner(mcmc_chain,bins=nbins,labels=[ p+" [\"]" for p in param_names], show_titles=True)
        plot.savefig(str(save_dir)+"/MCMC_multiplication_Df.png")
    title_plot = r"Combined $\Delta\phi$ from SDSSJ1433 lens mass model"
    if  KDE:
        from Plots.plotting_tools import plot_probability3D_KDE
        plot = plot_probability3D_KDE(Combined_PDF,Positions_KDE,labels=param_names,udm="arcsec^2",title=title_plot)
        plot.savefig(str(save_dir)+"/CombinedProbability_KDE.pdf") 
    else:
        from Plots.plotting_tools import plot_probability3D
        plot = plot_probability3D(Combined_PDF,Combined_bins,labels=param_names,udm="arcsec^2",title=title_plot)
        plot.savefig(str(save_dir)+"/CombinedProbability.pdf")

    print("Result directory:", str(save_dir))
    success(sys.argv[0])
