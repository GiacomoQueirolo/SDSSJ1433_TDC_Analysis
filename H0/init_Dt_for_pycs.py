# I want to have a standard way to obtain Dt
import sys,os
import numpy as np
import argparse as ap

from Utils.tools import present_program,success
from Plots.plotting_tools import get_median_with_error
from Utils.statistical_tools import marginalise_prob
from Utils.Dt_from_Df_reworked import Dt_XY
from Posterior_Analysis.tools_Post import default_cosmo
from H0.H0_Data import get_Df_post
from H0.H0_Combined_reworked import kwlnst

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Get prior value of time delay given dummy H0 and combined Df",
                               formatter_class=ap.RawTextHelpFormatter)
    
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of combined setting file)"
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default=".",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
    
    args      = parser.parse_args()
    df_name   = args.df_name
    combined_setting,PDF_Df,PDF_Df_bins = get_Df_post(df_name=df_name,link=False,overwrite=False,**kwlnst)
    Dfs_bin_densities =  marginalise_prob(PDF_Df,PDF_Df_bins)
    print("Using H0 "+str(default_cosmo.H0.value))
    dts = []
    for i in range(3):
        cnt = Dfs_bin_densities[i]
        bn  = PDF_Df_bins[i]
        Df_med,sig = get_median_with_error(cnt,bn,ret_str=False)
        dt_i =(-np.round(Dt_XY(Df_med,default_cosmo.H0.value,setting=combined_setting),2))
        print(dt_i)
        dts.append(dt_i)
    print(dts)
    success(sys.argv[0])
