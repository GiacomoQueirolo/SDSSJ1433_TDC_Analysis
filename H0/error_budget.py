import os,sys
import numpy as np
import argparse as ap

from H0.H0_Data import *

from Utils.tools import *
from Utils.get_res import load_whatever
from Utils.statistical_tools import estimate_sigma,estimate_median,quantiles,marginalise_prob

def get_error_budget(PH0_resdir,kwargs_dt,outprint=True,name_ph0_res="ph0_results",name="",**kwargs_df):
    sigma_dt = np.diag(kwargs_dt["cov"])**.5
    dts =kwargs_dt["mean"] 
    print(sigma_dt,dts,np.round(sigma_dt/dts,2))
    Dt_err_cntrb = np.mean(sigma_dt/dts)

    if outprint:
        str_i=f"relative error contribution for Dt: {np.round(Dt_err_cntrb*100,2)}%"
        if name:
            str_i = f"{name}: "+str_i
        print(str_i)
        
    PDF_Df,PDF_Df_bins = kwargs_df["PDF_Df"],kwargs_df["PDF_Df_bins"]

    dfs,sigma_df = [],[]
    bin_densities = marginalise_prob(PDF_Df,PDF_Df_bins)
    for i in range(3):
        bn  =  PDF_Df_bins[i]
        cnt = bin_densities[i]
        median = estimate_median([cnt],[bn])[0]
        sigmas = estimate_sigma([cnt],[bn],median=[median],averaged=True)[0]
        dfs.append(median)
        sigma_df.append(sigmas)
    dfs          = np.array(dfs)
    sigma_df     = np.array(sigma_df)
    print(sigma_df,dfs,np.round(sigma_df/dfs,2))
    Df_err_cntrb = np.mean(sigma_df/dfs)
    if outprint:
        str_i=f"relative error contribution for Df: {np.round(Df_err_cntrb*100,2)}%"
        if name:
            str_i = f"{name}: "+str_i
        print(str_i)
    
    PH0,H0 = load_whatever(f"{PH0_resdir}/{name_ph0_res}.data")
    h0_res,err_min,err_max = quantiles(PH0,H0,q=[0.16,.5,0.84],return_quantiles=False)
    
    sigmaH0  = np.mean([err_min,err_max])

    rel_err  = sigmaH0/h0_res
    print( [ np.round(np.sqrt(i**2+j**2),2) for i,j in zip(sigma_dt/dts,sigma_df/dfs)])

    print(sigmaH0,h0_res,np.round(rel_err,2))
    if outprint:
        str_i=f"Relative err. : {np.round(rel_err*100,2)}%"
        if name:
            str_i = f"{name}: "+str_i
        print(str_i)
        
    err_bdg =  {"Dt_err_cntrb":Dt_err_cntrb,"Df_err_cntrb":Df_err_cntrb,"rel_err":rel_err}

    with open(f"{PH0_resdir}/{name}error_budget.pkl","wb") as f:
        pickle.dump(err_bdg,f)

    print("The reason for which the rel. uncertainty of H0 is not = sqrt(rel.uncert.(dt)^2 +rel.uncert.(dphi)^2 )")
    print("is due to the fact that the 3 measurements are considered independent, and thus the uncertainty has a factor sqrt(3) of diff.")
    return err_bdg

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Combine the posterior for time delay and Fermat potential differences at the images position to constrain the Hubble parameter H0",
                               formatter_class=ap.RawTextHelpFormatter)
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of savedir for multiplied posteriors)"
    help_postdir   = "Name of the directory which will be containing the H0 combined posterior"
    parser.add_argument("-dPH0","--dir_PostH0",dest='dir_ph0', type=str,default="PH0",
                        metavar='dir_ph0', action='store',
                        help=help_postdir)
    parser.add_argument("-dtn","--Dt_name",dest='dt_name', type=str,default=".",
                        metavar='dt_name', action='store',
                        help=help_timedelay)
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default=".",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
    args      = parser.parse_args()
    dir_ph0   = args.dir_ph0
    dt_name   = args.dt_name
    df_name   = args.df_name
    res_dir   = "./results/"
    mkdir(res_dir)
    PH0_resdir = create_PH0resdir(res_dir=res_dir,dir_ph0=dir_ph0,verbose=False,overwrite=False)

    # default paths to directory for pycs, pycs config and lenstronomy
    kwpycs = {"pycs_path":"./time_J1433/","configpath":"/myconfig/"}
    kwlnst = {"lenstronomy_path":"./lens_J1433/"}

    # loading data
    config,Dt_res = get_Dt_post(dt_name=dt_name,PH0_resdir=PH0_resdir,overwrite=False,link=False,**kwpycs)
    combined_setting,PDF_Df,PDF_Df_bins = get_Df_post(df_name=df_name,PH0_resdir=PH0_resdir,overwrite=False,link=False,**kwlnst)
    #kwargs_df = {"combined_setting":combined_setting,"PDF_Df":PDF_Df,"PDF_Df_bins":PDF_Df_bins}
    # converting dt into kwargs
    kwargs_dt = get_kwdt(Dt_res)    
    print("All data collected")
    get_error_budget(PH0_resdir=PH0_resdir,kwargs_dt=kwargs_dt,PDF_Df=PDF_Df,PDF_Df_bins=PDF_Df_bins,outprint=True)
    success(sys.argv[0])
