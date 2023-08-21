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
import matplotlib
import numpy as np
import corner,pickle
import json,copy,time
import argparse as ap
import matplotlib.pyplot as plt

from H0_Combined_reworked import get_kwdt,get_PH0
from astropy.cosmology import Planck18,WMAP9,FlatLambdaCDM
from astropy.cosmology import FlatLambdaCDM
from lens_J1433.tools import *
from lens_J1433.plot_H0 import plot_H0
from tools_H0 import H0_Res
from lens_J1433.Dt_from_Df_reworked import *
from lens_J1433.statistical_tools import quantiles
from H0_Data import *

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Test the effect of cosmology in k for H0 ",
                               formatter_class=ap.RawTextHelpFormatter)
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of savedir for multiplied posteriors)"
    help_postdir   = "Name of the directory which will be containing the H0 combined posterior"
    help_overwrite = "Overwrite previous result. If False and I find the same result directory, I will create a new one with today's date."
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
    parser.add_argument('-ovw','--overwrite',help=help_overwrite,
                        dest="overwrite", 
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
    verbose   = args.verbose
    overwrite = args.overwrite

    res_dir   = "./results_testK/"
    mkdir(res_dir)
    PH0_resdir = create_PH0resdir(res_dir=res_dir,dir_ph0=dir_ph0,verbose=verbose,overwrite=overwrite)

    # loading data
    Dt_res = get_Dt_post(dt_name=dt_name,PH0_resdir=PH0_resdir,overwrite=overwrite,**kwpycs)
    combined_setting,PDF_Df,PDF_Df_bins = get_Df_post(df_name=df_name,PH0_resdir=PH0_resdir,overwrite=overwrite,**kwlnst)
    write_readme(path=PH0_resdir,dt_name=dt_name,df_name=df_name,**kwpycs,**kwlnst)
    # converting dt into kwargs
    kwargs_dt = get_kwdt(Dt_res)    

    print("All data collected")

    # Computing H0 posterior    
    H0_sampled = np.arange(h0min,h0max,h0step)
    
    color = ["k","b","r"]
    EdS_cosmo = FlatLambdaCDM(H0=100,Om0=1,name="EdS")
    f, ax = plt.subplots(1, 1, figsize=(18, 10))
    for ic,cosmo in enumerate([Planck18,WMAP9,EdS_cosmo]):
        PH0,H0 = get_PH0(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,setting=combined_setting,cosmo=cosmo)
        h0_res,err_min,err_max = quantiles(PH0,H0,q=[0.16,.5,0.84],return_quantiles=False)
        H0_res = H0_Res(h0_res,[err_min,err_max])
        print("Cosmo: ", cosmo.name)
        print("analytical: ", H0_res)
        ax = plot_H0(H0,PH0,add_mode=False,color=color[ic],return_plot=True,ax=ax)
    figname =f"{PH0_resdir}/PH0_dtf_cosmo.pdf"
    plt.savefig(figname)
    plt.close()
    print(f"Saved {figname}")
    
    Omega_m_smpl = np.linspace(0,1,100)
    PH0_Om  = []
    for Om in Omega_m_smpl:
        cosmo  = FlatLambdaCDM(H0=Planck18.H0,Om0=Om,name=f"Omega_m:{np.round(Om,2)}")
        PH0,H0 = get_PH0(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,setting=combined_setting,cosmo=cosmo)
        PH0_Om.append(PH0)
    Om,Ho = np.meshgrid(Omega_m_smpl,H0)


    plt.pcolormesh(Om,Ho,np.transpose(PH0_Om))
    savename = f"{PH0_resdir}/PH0_OmegaM.pdf"
    plt.savefig(savename)
    print(f"Saved {savename}")
    PH0_OM = [Om,Ho,np.transpose(PH0_Om)]
    with open(f"{PH0_resdir}/PH0_OM.pkl","wb") as f:
        pickle.dump(PH0_OM,f)
    success(sys.argv[0])

    
