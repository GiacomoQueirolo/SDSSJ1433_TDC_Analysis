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
from datetime import datetime
from scipy.stats import multivariate_normal

from Utils.tools import *
from Plots.plot_H0 import plot_H0
from H0.tools_H0 import H0_Res
from Utils.Dt_from_Df_reworked import *
from Utils.statistical_tools import quantiles
from Utils.combinedsetting_class import combined_setting


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

def get_PH0(Dt_kw,Df_dens,Df_dens_bins,setting,H0=np.arange(50,100,.1)):
    PH0     = []
    k       = k_analytical(setting) 
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



def write_readme(path,dt_path,df_path):
    # save a txt file in the final path specifing the path
    # used for to find the dt and df data for the posterior
    with open(str(path)+"/paths_to_data_used.txt","w") as f:
        f.write("Time delay posterior used:\n")
        f.write(str(dt_path)+"\n")
        f.write("Fermat potential posterior used:\n")
        f.write(str(df_path)+"\n")
    return 0

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Combine the posterior for time delay and Fermat potential differences at the images position to constrain the Hubble parameter H0",
                               formatter_class=ap.RawTextHelpFormatter)
    help_timedelay = "Name of the directory containing the posterior of the difference of Time Delay"
    help_fermatpot = "Name of the directory containing the posterior of the difference of Fermat Potential"
    help_postdir   = "Name of the directory which will be containing the H0 combined posterior"
    help_overwrite = "Overwrite previous result. If False and I find the same result directory, I will create a new one with today's date."
    parser.add_argument("-dPH0","--dir_PostH0",dest='dir_ph0', type=str,default="PH0",
                        metavar='dir_ph0', action='store',
                        help=help_postdir)
    parser.add_argument("-ddt","--dir_Dtimedelay",dest='dir_dt', type=str,default=".",
                        metavar='dir_td', action='store',
                        help=help_timedelay)
    parser.add_argument("-ddf","-dir_Dfermat",dest='dir_df', type=str,default=".",
                        metavar='dir_df', action='store',
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
    args    = parser.parse_args()
    dir_ph0 = args.dir_ph0
    dir_dt  = args.dir_dt
    dir_df  = args.dir_df
    h0max   = args.h0max
    h0min   = args.h0min
    h0step  = args.h0step
    verbose = args.verbose
    overwrite = args.overwrite
    
    res_dir = "./results/"
    mkdir(res_dir)
    PH0_resdir = res_dir+"/"+dir_ph0
    if os.path.isfile(PH0_resdir):
        if os.listdir(PH0_resdir):
            if verbose:
                print(f"Previous results found in {PH0_resdir}")
            if overwrite:
                if verbose:
                    print(f"Overwriting them (actually moving them to old_{PH0_resdir})")
                os.rename(PH0_resdir, "old_"+PH0_resdir)
            else:
                td = datetime.today().strftime("%y%m%d")
                if verbose:
                    print(f"Changing result directory name: {PH0_resdir}->{PH0_resdir}_{td}")
                PH0_resdir = PH0_resdir+"_"+td
    mkdir(PH0_resdir)
    print("Saving results in: ",PH0_resdir)
    dt_resdir  = PH0_resdir+"/Dt_post"
    os.symlink(os.getcwd()+"/"+str(dir_dt),dt_resdir)
    print("For time delay using posterior : ",str(dir_dt))
    df_resdir  = PH0_resdir+"/Df_post"
    os.symlink(os.getcwd()+"/"+str(dir_df),df_resdir)
    print("For fermat pot. using posterior: ",str(dir_df))

    with open(df_resdir+"/Combined_PDF.pkl","rb") as f:
        PDF_Df = np.array(pickle.load(f))
    with open(df_resdir+"/Combined_PDF_bins.pkl","rb") as f:
        PDF_Df_bins = np.array(pickle.load(f))
        
    with open(dt_resdir+"/marginalisation_spline_sigma_0.50_combined.pkl","rb") as f:
        Dt_res = pickle.load(f)
        
    kwargs_dt = get_kwdt(Dt_res)

    print("All data collected")

    ### the combined setting is only in the fermat pot dir
    with open(df_resdir+"/combined_setting.pkl","rb") as f:
        combined_setting = pickle.load(f)
    
    H0_sampled = np.arange(h0min,h0max,h0step)
    if verbose:
        begin = time.time()
    PH0,H0 = get_PH0(Dt_kw=kwargs_dt,Df_dens=PDF_Df,Df_dens_bins=PDF_Df_bins,H0=H0_sampled,setting=combined_setting)
    if verbose:
        print("Time passed: ",time.time()-begin)
    
    with open(PH0_resdir+"/ph0_results.data","wb") as f:
        pickle.dump([PH0,H0],f)
    
    h0_res,err_min,err_max = quantiles(PH0,H0,q=[0.16,.5,0.84],return_quantiles=False)
    H0_res = H0_Res(h0_res,[err_min,err_max])
    print("analytical: ", H0_res)
    plot_H0(H0,PH0,figname=PH0_resdir+"/PH0_dtf.pdf",add_mode=False)
    write_readme(PH0_resdir,dt_resdir,df_resdir)
    success(sys.argv[0])


"""

def get_normalisation(in_bins,Dt_kw,H0_range):
    bins = np.array(copy.deepcopy(in_bins))
    dH0 = H0_range[1]-H0_range[0] #is this necessary?
    kwargs_df = {"mean":Df_XY(copy.deepcopy(Dt_kw["mean"]),H0_range[0]),
     "cov":cov_Df(copy.deepcopy(Dt_kw["cov"]),H0_range[0])}
    norm = get_analytic_density(bins=bins,**kwargs_df)
    for i in range(1,len(H0_range)):
        kwargs_df = {"mean":Df_XY(copy.deepcopy(Dt_kw["mean"]),H0_range[i]),
             "cov":cov_Df(copy.deepcopy(Dt_kw["cov"]),H0_range[i])}
        norm+= get_analytic_density(bins=bins,**kwargs_df)*dH0
    print("normalisation done")
    return norm

def get_bin_vol(in_bins):
    bins =copy.deepcopy(in_bins)
    dim = np.shape(bins)[0]
    print(dim)
    bin_vol = 1
    for d in range(dim):
        bin_lenght=bins[d][1]-bins[d][0] 
        for n in range(5):
            rnd_index = np.random.randint(2,len(bins[d]))
            test_lenght =bins[d][rnd_index]-bins[d][rnd_index-1]
            if np.abs(bin_lenght-test_lenght)>1e-7:
                raise RuntimeError("Not all bins have the same size")
        bin_vol *= bin_lenght
    # we get a "volume" (vol in 3D, area in 2D, lenght in 1D) 
    # assuming that all bins are the same
    return bin_vol

def get_PH0_plot(H0=np.arange(50,100,.1),
            Dt_kw = kwargs_dt,
            #Df    = mcmc_df,
            Dens_f = PDF_Df,
            nd_bins = PDF_Df_bins,
            savefigs=False):
    
    PH0 = []
    Dens_f_trsf_norm = get_normalisation(in_bins=copy.deepcopy(nd_bins),\
                                             Dt_kw = copy.deepcopy(Dt_kw),\
                                             H0_range=H0)
    for h0 in H0:
        kwargs_df = {"mean":Df_XY(copy.deepcopy(Dt_kw["mean"]),h0),
             "cov":cov_Df(copy.deepcopy(Dt_kw["cov"]),h0)}
        
        
        Dens_f_trsf = get_analytic_density(bins=copy.deepcopy(nd_bins),\
                                           norm=Dens_f_trsf_norm,\
                                           **kwargs_df)
        
        Dens_tot = copy.deepcopy(Dens_f)*Dens_f_trsf
        # the integration in this case is nothing else then the sum over every bin, ie:
        Darea_bins = get_bin_vol(copy.deepcopy(nd_bins))
        P_h0 = np.sum(Dens_tot)*Darea_bins
        if np.isnan(P_h0):
            P_h0 = 0.
        PH0.append(P_h0)
        
    if np.sum(PH0)!=0:
        print("Evidence: ",np.sum(PH0))
        PH0 = PH0/np.sum(PH0)
    else:
        print("WARNING: sum of PH0 == 0")
    return np.array(PH0),H0

from importlib import reload 
import plot_H0
plot_H0 = reload(plot_H0)#.plot_H0
plot_H0 = plot_H0.plot_H0
#from Plots.plot_H0 import plot_H0
plot_H0(H0,PH0,figname=res_dir+"/PH0_dtf"+CP+".png",add_mode=False,transparent=False)

"""

