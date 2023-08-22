#!/usr/bin/env python
# coding: utf-8


"""
Inverse of Dft_combo :

The idea is to convert the 3D Dt distribution
using some given test values of H0 into Df distr.
and see how well it agrees with the 3D df distribution

""";


# $$ \large P (H_0)\ =\ \int d\Delta \phi P^{\text{lens model}}(\Delta \phi)* \frac{P^{\text{lightcurves}}(\Delta \phi (\Delta t)|H_0)}{\int dH_0 ' P^{\text{lightcurves}}(\Delta \phi (\Delta t)|H_0') }*Prior(H_0) $$

# ### What we want to obtain is:
# 
# $$ \large P (H_0)\ =\ \int\int\int d\Delta \phi_{AB}d\Delta \phi_{AC}d\Delta \phi_{BC} P_{\Delta \phi\  meas.}(\Delta \phi_{AB},\Delta \phi_{AC},\Delta \phi_{BC})* \frac{P_{\Delta t\ transf.}(\Delta \phi_{AB},\Delta \phi_{AC},\Delta \phi_{BC}|H_0)*Prior(H_0)}{\int dH_0 ' P_{\Delta t \  transf.}(\Delta \phi_{AB},\Delta \phi_{AC},\Delta \phi_{BC} | H_0') } $$

import os,sys
import numpy as np
import pickle
import copy,time
import argparse as ap
import multiprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from astropy.cosmology import Planck18
from astropy.cosmology import FlatLambdaCDM
from scipy.stats import multivariate_normal

from Utils.tools import *
from Plots.plot_H0 import plot_H0
from H0.tools_H0 import H0_Res
from Utils.Dt_from_Df_reworked import *
from Utils.statistical_tools import quantiles
from Utils.combinedsetting_class import combined_setting
from H0_Data import get_Df_post,get_Dt_post,create_PH0resdir,write_readme

def get_kwdt(Dt_res):
    #MOD_MULTIVAR
    # we want to consider the error correlation, see notes 12th Nov.
    cov_rnd   = np.cov(Dt_res.err_distr)
    cov_sys   = (np.array(Dt_res.error.sys)**2)*np.identity(len(Dt_res.error.sys))
    cov_dt    = cov_sys + cov_rnd
    kwargs_dt = {"mean":Dt_res.results,"cov":cov_dt}
    return kwargs_dt

def get_bin_center(bins):
    bins_centers = []
    for i in range(len(bins)):
        bins_centers_i = []
        for j in range(len(bins[i])-1):
            bins_centers_i.append((bins[i][j]+bins[i][j+1])/2.)
        bins_centers.append(bins_centers_i)
    return  np.array(bins_centers)

def get_analytic_density(mean,cov,bins,norm=1):
    # assuming gaussian probability density function
    center_bins    = get_bin_center(bins) 
    grid_of_points = np.transpose(np.meshgrid(*center_bins))
    density        = multivariate_normal.pdf(grid_of_points,mean,cov)*norm
    return density

def get_PH0(Dt_kw,Df_dens,Df_dens_bins,setting,H0=np.arange(50,100,.1),cosmo=Planck18):
    PH0     = []
    k       = k_analytical(setting,cosmo=cosmo) 
    for h0 in H0:
        # copy the fix elements in order to be sure 
        # not to modify them during the computation
        Df_dens_lens = copy.deepcopy(Df_dens)
        bins_Df      = copy.deepcopy(Df_dens_bins)
        dt_mean      = copy.deepcopy(Dt_kw["mean"])
        dt_cov       = copy.deepcopy(Dt_kw["cov"])
        kw_df_trnsf  = {"mean":Df_XY(Dt_XY=dt_mean,H0=h0,k=k),
                        "cov":cov_Df(cov_Dt=dt_cov,H0=h0,k=k)}
        Df_dens_trsf = get_analytic_density(bins=bins_Df,**kw_df_trnsf)
        Df_dens_comb = Df_dens_lens*Df_dens_trsf
        # the integration in this case is nothing else then the sum over every bin, ie:
        P_h0 = np.sum(Df_dens_comb)
        if np.isnan(P_h0):
            print("Warning, nan P_h0")
            P_h0 = 0.
        PH0.append(P_h0)
    if np.sum(PH0)!=0:
        PH0 = PH0/np.sum(PH0)
    else:
        print("WARNING: sum of PH0 == 0")
    return np.array(PH0),H0

def _pll_PH0i_OmegaM(args):
    return _PH0i_OmegaM(*args)

def _PH0i_OmegaM(Omi,Dt_kw,Df_dens,Df_dens_bins,setting,H0,cosmo0):
    cosmoi  = FlatLambdaCDM(Om0=Omi,H0=cosmo0.H0,Tcmb0=cosmo0.Tcmb0)
    PH0i,H0 = get_PH0(Dt_kw=Dt_kw,Df_dens=Df_dens,Df_dens_bins=Df_dens_bins,setting=setting,H0=H0,cosmo=cosmoi)
    #PH0i,H0 = get_PH0_nosett(Dt_kw=Dt_kw,Df_dens=Df_dens,Df_dens_bins=Df_dens_bins,\
    #                           z_lens=setting.z_lens,z_source=setting.z_source,H0=H0,cosmo=cosmoi)
    return PH0i,H0
    
def get_PH0_marg_OmegaM(Dt_kw,Df_dens,Df_dens_bins,setting,H0=np.arange(50,100,.1),cosmo0 = Planck18,Omega_M=np.arange(0,1.002,.002),parall=False):
    if parall:
        def _get_PH0_marg_OmegaM(Omi):
            cosmoi  = FlatLambdaCDM(Om0=Omi,H0=cosmo0.H0,Tcmb0=cosmo0.Tcmb0)
            PH0i,H0i = get_PH0(Dt_kw=Dt_kw,Df_dens=Df_dens,Df_dens_bins=Df_dens_bins,setting=setting,H0=H0,cosmo=cosmoi)
            return PH0i,H0i
        with multiprocess.Pool() as pool:
            #PH0s,verify_H0s = [pool.starmap(_pll_PH0i_OmegaM,(Omi,Dt_kw,Df_dens,Df_dens_bins,setting,H0,cosmo0)) for Omi in Omega_M]
            PH0s,verify_H0s = zip(*pool.map(_get_PH0_marg_OmegaM,Omega_M))
    else:    
        PH0s = []
        verify_H0s = [H0]
        for Omi in Omega_M:
            PH0i,H0 = _PH0i_OmegaM(Omi,Dt_kw,Df_dens,Df_dens_bins,setting,H0,cosmo0)
            PH0s.append(PH0i)
            verify_H0s.append(H0)
    for H0i in verify_H0s:
        if np.array(H0i).tolist()!=np.array(verify_H0s[0]).tolist():
            raise RuntimeError("Something went wrong with the sampling")
    H0 = verify_H0s[0]
    # PH0_marg = integ ( d Omega_M * P(H0,Omega_m)) = sum (d Omega_M[const] * P(H0,Omega_m)) = dOmega_M*sum(P(H0,Omega_m))
    Marg_PH0 = np.sum(PH0s,axis=0)*0.001
    norm_PH0 = Marg_PH0/np.sum(Marg_PH0)
    return norm_PH0,H0,PH0s

# default paths to directory for pycs, pycs config and lenstronomy
kwpycs = {"pycs_path":"./time_J1433/","configpath":"/myconfig/"}
kwlnst = {"lenstronomy_path":"./lens_J1433/"}


if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Combine the posterior for time delay and Fermat potential differences at the images position to constrain the Hubble parameter H0",
                               formatter_class=ap.RawTextHelpFormatter)
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"
    #help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of savedir for multiplied posteriors)"
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of combined setting file)"
    help_postdir   = "Name of the directory which will be containing the H0 combined posterior"
    help_overwrite = "Overwrite previous result. If False and I find the same result directory, I will create a new one with today's date."
    help_margOm    = "Consider marginalised Posterior with respect to Omega_m"
    parser.add_argument("-dPH0","--dir_PostH0",dest='dir_ph0', type=str,default="PH0",
                        metavar='dir_ph0', action='store',
                        help=help_postdir)
    parser.add_argument("-dtn","--Dt_name",dest='dt_name', type=str,default=".",
                        metavar='dt_name', action='store',
                        help=help_timedelay)
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default=".",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
    parser.add_argument("-h0min",type=int, dest="h0min", default=50,
                    help="Minumum value for H0 sampling")
    parser.add_argument("-h0max",type=int, dest="h0max", default=100,
                    help="Maximum value for H0 sampling")
    parser.add_argument("-h0step",type=float, dest="h0step", default=.1,
                    help="Step for H0 sampling")
    parser.add_argument('-mOm','--marg_Om',help=help_margOm,
                        dest="marg_Om", 
                        default=False,action="store_true")
    parser.add_argument('-ovw','--overwrite',help=help_overwrite,
                        dest="overwrite", 
                        default=False,action="store_true")
    parser.add_argument('-pll','--parallel',help="Use parallel computing",
                        dest="parellel", 
                        default=False,action="store_true")
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args      = parser.parse_args()
    dir_ph0   = args.dir_ph0
    dt_name   = args.dt_name
    df_name   = args.df_name
    h0max     = args.h0max
    h0min     = args.h0min
    h0step    = args.h0step
    marg_Om   = args.marg_Om 
    verbose   = args.verbose
    parallel  = args.parellel
    overwrite = args.overwrite

    
    res_dir   = "./results/"
    mkdir(res_dir)
    PH0_resdir = create_PH0resdir(res_dir=res_dir,dir_ph0=dir_ph0,verbose=verbose,overwrite=overwrite)

    # loading data
    Dt_res = get_Dt_post(dt_name=dt_name,PH0_resdir=PH0_resdir,overwrite=overwrite,**kwpycs)
    combined_setting,PDF_Df,PDF_Df_bins = get_Df_post(df_name=df_name,PH0_resdir=PH0_resdir,overwrite=overwrite,**kwlnst)
    write_readme(path=PH0_resdir,dt_name=dt_name,df_name=df_name,log_command=True,**kwpycs,**kwlnst)
    # converting dt into kwargs
    kwargs_dt = get_kwdt(Dt_res)    

    print("All data collected")

    # Computing H0 posterior    
    H0_sampled = np.arange(h0min,h0max,h0step)
    if verbose:
        begin = time.time()
    PH0,H0 = get_PH0(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,setting=combined_setting,cosmo=Planck18)
    if verbose:
        print("Time passed: ",time.time()-begin)
    
    # saving result
    with open(PH0_resdir+"/ph0_results.data","wb") as f:
        pickle.dump([PH0,H0],f)
    # plotting and outputting result
    h0_res,err_min,err_max = quantiles(PH0,H0,q=[0.16,.5,0.84],return_quantiles=False)
    H0_res = H0_Res(h0_res,[err_min,err_max])
    print("analytical: ", H0_res)
    plot_H0(H0,PH0,figname=PH0_resdir+"/PH0_dtf.pdf",add_mode=False)
    
    if marg_Om:
        f, ax = plt.subplots(1, 1, figsize=(18, 10))
        ax = plot_H0(H0,PH0,add_mode=False,color="blue",return_plot=True,ax=ax)
        handles = [Patch(facecolor="blue",label=r"$\Omega_{m,0}$"+f"={np.round(Planck18.Om0,3)}")]
        # uniform Om = U[0,1]
        marg_PH0,H0,PH0_2D = get_PH0_marg_OmegaM(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,setting=combined_setting,cosmo0=Planck18,parall=parallel)
        with open(f"{PH0_resdir}/marg_ph0_results.data","wb") as f:
            pickle.dump([marg_PH0,H0],f)    
        ax = plot_H0(H0,marg_PH0,ax=ax,return_plot=True,color="cyan",add_mode=False)
        handles.append(Patch(facecolor="cyan",label=r"$\Omega_{m,0}$= $\mathcal{U}$(0.,1.)"))
        # the followings are inspired by TDCOSMO XIII, Section 6
        # uniform Om = U[0.05,.5]
        marg_PH0,H0,PH0_2D = get_PH0_marg_OmegaM(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,\
                            Omega_M=np.arange(0.05,0.502,0.001),setting=combined_setting,cosmo0=Planck18,parall=parallel)
        with open(f"{PH0_resdir}/marg_ph0_results_OmU005_05.data","wb") as f:
            pickle.dump([marg_PH0,H0],f)    
        ax = plot_H0(H0,marg_PH0,ax=ax,return_plot=True,color="darkorange",add_mode=False)
        handles.append(Patch(facecolor="darkorange",label=r"$\Omega_{m,0}$= $\mathcal{U}$(0.05,0.5)"))
        # uniform Om = N[0.334, 0.018]
        marg_PH0,H0,PH0_2D = get_PH0_marg_OmegaM(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,\
                            Omega_M=np.random.normal(0.334,0.018,1000),setting=combined_setting,cosmo0=Planck18,parall=parallel)
        with open(f"{PH0_resdir}/marg_ph0_results_OmN0334_0018.data","wb") as f:
            pickle.dump([marg_PH0,H0],f)    
        ax = plot_H0(H0,marg_PH0,ax=ax,return_plot=True,color="yellow",add_mode=False)
        handles.append(Patch(facecolor="yellow",label=r"$\Omega_{m,0}$= $\mathcal{U}$(0.05,0.5)"))
        
        ax.legend(handles=handles)
        figname =f"{PH0_resdir}/PH0_dtf_cosmo_marg_extended.pdf"
        plt.savefig(figname)
        plt.close()
        print(f"Saved {figname}")
    
        
    success(sys.argv[0])