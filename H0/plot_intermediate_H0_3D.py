#!/usr/bin/env python

"""
Plan to obtain intermediate H0 from single df of single filter model of lenstronomy and overplot the posterior
""";


import os,sys,copy
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm as normal
from astropy.cosmology import Planck18

from tools_H0 import *
from tools import *
from plot_H0 import plot_H0
from get_res import get_combined_Df
from plotting_tools import base_colors
from fermat_pot_analysis import get_mcmc_Df
from combinedsetting_class import combined_setting
from statistical_tools import quantiles,marginalise_prob
from Dt_from_Df_reworked import k_analytical,Df_XY,cov_Df

from H0_Data import get_Dt_post,create_PH0resdir
from H0_Combined_reworked import get_PH0,kwpycs,kwlnst,get_kwdt

lnstr_path       = f'./{kwlnst["lenstronomy_path"]}/'

def get_bin_center(bins):
    if len(np.shape(bins))==1: #so it works for 1D too
        bins=[bins]
    bins_centers = []
    for i in range(len(bins)):
        bins_centers_i = []
        for j in range(len(bins[i])-1):
            bins_centers_i.append((bins[i][j]+bins[i][j+1])/2.)
        bins_centers.append(bins_centers_i)
    return np.array(bins_centers)

def get_analytic_density_1D(mean,err,bins,norm=None):
    center_bins = get_bin_center(bins) 
    #grid_of_points = np.transpose(np.meshgrid(*center_bins))
    dens  = normal.pdf(center_bins,loc=mean,scale=err)
    if type(norm) is float or type(norm) is int:
        dens/=norm
    elif type(norm) is bool:
        if norm:
            dens/=np.sum(dens)
    return dens

def err_Df(err_Dt,i_dim,H0,k):
    # Convert sigma dt in sigma dphi 
    #-> simply the sqrt of the covariance
    cov_Df_ = cov_Df(err_Dt**2,H0,k)
    err_Df_ = np.sqrt(cov_Df_)
    return err_Df_[i_dim][i_dim]

def get_PH0_1D(Dt_kw,Df_dens,Df_dens_bins,setting,i_dim,H0=np.arange(50,100,.1),cosmo=Planck18):
    PH0     = []
    k       = k_analytical(setting,cosmo=cosmo) 
    for h0 in H0:
        # copy the fix elements in order to be sure 
        # not to modify them during the computation
        Df_dens_lens = copy.deepcopy(Df_dens)
        bins_Df      = copy.deepcopy(Df_dens_bins)
        dt_mean      = copy.deepcopy(Dt_kw["mean"])
        dt_cov       = copy.deepcopy(Dt_kw["cov"])
        kw_df_trnsf  = {"mean" : Df_XY(Dt_XY=dt_mean,H0=h0,k=k)[i_dim],
                        "err"  : err_Df(err_Dt=dt_cov,i_dim=i_dim,H0=h0,k=k)}
        Df_dens_trsf = get_analytic_density_1D(bins=bins_Df,**kw_df_trnsf)
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


if __name__=="__main__":
    ################################
    name_prog = sys.argv[0]
    present_program(name_prog)
    ################################
    parser = argparse.ArgumentParser(description="Temporary plot of the posterior of H0 from the given lens models divided in lcs couples (AB,AC,BC-AD)")
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of combined setting file)"
    help_numbins   = "Number of bins per dimension for the histog. sampling of the Df (be careful with it! too many bins can be catastrophic)"
    help_dir       = "In which PH0 directory to save it (should be the same as the one used for H0_Combined_reworked)"
    parser.add_argument("-dtn","--Dt_name",dest='dt_name', type=str,default=".",
                        metavar='dt_name', action='store',
                        help=help_timedelay)
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                    help=help_numbins)
    parser.add_argument("-dir","--directory",type=str, dest="dir", default="PH0",
                    help=help_dir)
    parser.add_argument("-sl","--simple_legend",dest="simple_legend", default=False,action="store_true",
                        help="Draw a simplified legend with only the name of the filters")
    parser.add_argument("-AD", action="store_true", dest="AD",default=False,
                        help="Consider AD couple instead of BC")
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")

    parser.add_argument("-cs","-CombSett", action="store_true", dest="CombSett",default=False,
                        help="Used combined posterior of the Fermat Potential between the settings")
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default="",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    # bins for the fermat pot
    bins = args.nbins 
    dir  = args.dir
    verbose         = args.verbose
    simple_legend   = args.simple_legend
    setting_names   = args.SETTING_FILES
    BC              = not bool(args.AD)
    CombSett        = args.CombSett
    df_name         = args.df_name
 
    PH0_resdir = create_PH0resdir(dir_ph0=dir,verbose=verbose,overwrite=False)

    if CombSett and df_name=="":
        filters       = [get_filter(st) for st in setting_names]
        backup_path  = "backup_results"
        main_savedir = "PDF_multiplication_ABC"
        if not BC:
            main_savedir+="D"
        save_dir      = create_dir_name(setting_names,save_dir=main_savedir,dir_name=lnstr_path,backup_path=backup_path,copy_settings=True)
        cmb_sett_name = combined_setting("",1,0,filters,setting_names,savedir=save_dir).gen_cmb_sett_name(cmb_sett_name=None)
        try:
            combined_setting = get_combined_setting_module(combined_sett,main_dir=lnstr_path)
            print(f"Recovered combined setting from setting list: {cmb_sett_name}")
        except:
            print("If you want to use the combined setting, either 'df_name' have to be given or the setting names must be in the right order") 
            exit()
        df_name = cmb_sett_name
    if not CombSett and df_name!="":
        raise ValueError("'CombSett' is False, but the 'df_name' is defined, something is wrong")
        
    dt_name      = args.dt_name
    Dt_res       = get_Dt_post(dt_name=dt_name,link=False,**kwpycs)
    kwargs_dt    = get_kwdt(Dt_res)    
    savefig_path = f"{lnstr_path}/backup_results/Post_H0/{dir}"
    mkdir(savefig_path)
         
    legend_elements = []
    labels_lcs = ["AB","AC","BC"]
    if not BC:
        labels_lcs[-1] = "AD"
    f,axes = plt.subplots(2,2,figsize=(12,12)) 

    if CombSett:
        combined_setting, Combined_PDF,Combined_bins = get_combined_Df(combined_sett=df_name,main_dir=lnstr_path)
        #Note: doesn't work with KDE
        marg_Combined_PDF = marginalise_prob(Combined_PDF,Combined_bins)
        print(f"For Fermat Potential using posterior: {combined_setting.cmb_sett_name}")

    for dim_i in range(3): # each dimension/quadrant
        if dim_i==0:
            ax = axes[0][0]
        elif dim_i==1:
            ax = axes[1][0]
        elif dim_i==2:
            ax = axes[1][1]
        else:
            raise RuntimeError("")
        if not CombSett:
            for i,sets in enumerate(setting_names): #each setting
                print_setting(sets)
                # a bit hacky way of getting the setting module, but it works
                orig_wd = os.getcwd()
                os.chdir(lnstr_path)
                sett  = get_setting_module(sets,1)
                Dfi   = get_mcmc_Df(sett,noD=BC)
                os.chdir(orig_wd)
                
                mcmc_Df = np.transpose(Dfi) #shape: 3, len(mcmc)
                if not BC:
                    print("WARNING: Given the low S/N of image D, I will discard here and instead consider Delta BC")    
                Post_Df,Post_Df_bins = np.histogramdd(mcmc_Df[dim_i].T,bins=bins,density=True)
                PH0,H0 = get_PH0_1D(Df_dens=Post_Df,Df_dens_bins=Post_Df_bins,Dt_kw=kwargs_dt,setting=sett,i_dim=dim_i)  
                ax = plot_H0(H0,PH0,add_mode=False,color=base_colors[i],return_plot=True,ax=ax,title=labels_lcs[dim_i])
                if dim_i==1:
                    for lg in ax.get_legend().legendHandles:
                        legend_elements.append(lg)
                ax.get_legend().remove()
                legend_elements.append(strip_setting_name(sett,filter=simple_legend))
        else:
            PH0,H0 = get_PH0_1D(Df_dens=marg_Combined_PDF[dim_i],Df_dens_bins=Combined_bins[dim_i],i_dim=dim_i,Dt_kw=kwargs_dt,setting=combined_setting)
            ax = plot_H0(H0,PH0,add_mode=False,return_plot=True,ax=ax,title=labels_lcs[dim_i])
            if dim_i==1:
                for lg in ax.get_legend().legendHandles:
                    legend_elements.append(lg)
            ax.get_legend().remove()
    axdel=axes[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    f.suptitle("Compare $H_0$ posterior from individual filter's lens models for different LCs combinations")
    plt.tight_layout()
    output_name = f"{PH0_resdir}/compare_H0_3D.pdf"
    if CombSett:
        output_name = output_name.replace(".pdf","_CombSett.pdf")
    if not BC:
        output_name = output_name.replace(".pdf","_ABCD.pdf")
    f.savefig(output_name)
    print(f"Created {output_name}")
    success(name_prog)    

