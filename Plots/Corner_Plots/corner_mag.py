# corner plot fermat pot and mag ratio
import json
import argparse
import numpy as np
from copy import copy
from corner import corner
import matplotlib.pyplot as plt

from Utils.tools import *
from Utils.get_res import *
from Data.conversion import qphi_from_e1e2
from Posterior_analysis.mag_remastered import labels as  labels_mag

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Plot the corner plot of mag ratios and all other lens params")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    settings=args.SETTING_FILES
    for sett in settings:
        mcmc_path = get_savemcmcpath(sett)
        sett_str  = strip_setting_name(sett)
        mag = np.array(load_whatever(mcmc_path+"/mcmc_mag_rt_"+sett_str+".json"))
        # mag is already a ratio
        smpl   = get_mcmc_smpl(sett)
        prm    = get_mcmc_prm(sett)
        e1     = smpl.T[prm.index("e1_lens0")]
        e2     = smpl.T[prm.index("e2_lens0")]
        q,phi  = qphi_from_e1e2(e1,e2,ret_deg=True)
        #q,phi  = np.array(q).reshape(len(q),1),np.array(phi).reshape(len(phi),1)
        smplT  = copy(smpl.T)
        smplT[prm.index("e1_lens0")] = q
        smplT[prm.index("e2_lens0")] = phi
        smpl   = np.hstack([mag,smplT.T])
        prm[prm.index("e1_lens0")] = "q_lens0"
        prm[prm.index("e2_lens0")] = "phi_lens0"        
        # truths from   "from_dmag_to_FR.py" in pycs with results from J1433_testFR
        print("truths from 'from_dmag_to_FR.py' in pycs with results from J1433_testFR")
        truths = [1.26,-0.85,-0.675 ,*np.median(smpl.T[3:],1).tolist()]
        labels = [*labels_mag,*prm]
        plot = corner(smpl,truths=truths, labels=labels,show_titles=True)
        plot.savefig(get_savefigpath(sett)+"/corner_mag.pdf")
        print("Saved plot "+get_savefigpath(sett)+"/corner_mag.pdf")
    success(sys.argv[0])

