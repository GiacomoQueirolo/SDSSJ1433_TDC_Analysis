# corner plot fermat pot and mag ratio
import json
import argparse
import numpy as np
from corner import corner
import matplotlib.pyplot as plt

from Utils.tools import *
from Utils.get_res import *
from Posterior_analysis.fermat_pot_analysis import labels_Df
from Posterior_analysis.mag_remastered import labels as  labels_mag

            
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Plot the corner plot of the Df and mag ratios")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    
    settings=args.SETTING_FILES
    for sett in settings:
        mcmc_path = get_savemcmcpath(sett)
        sett_str  = strip_setting_name(sett)
        mag    = np.array(load_whatever(mcmc_path+"/mcmc_mag_rt_"+sett_str+".json"))
        # mag is already a ratio
        fermat = np.array(load_whatever(mcmc_path+"/mcmc_ordered_fermat_"+sett_str+".json"))
        # fermat is still to be taken as a difference
        fermat = fermat.T[1:] - fermat.T[0]
        fermat = fermat.T
        mag_fermat = np.hstack((mag,fermat))
        plot = corner(mag_fermat,labels=[*labels_mag,*labels_Df],show_titles=True)
        plot.savefig(get_savefigpath(sett)+"/corner_mag_df.png")
        print("Saved plot "+get_savefigpath(sett)+"/corner_mag_df.png")
    success(sys.argv[0])

