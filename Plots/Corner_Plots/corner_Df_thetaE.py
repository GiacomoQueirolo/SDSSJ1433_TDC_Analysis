#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import importlib
import os,sys
import json
import argparse
import corner
from copy import copy

from Utils.tools import *
from Utils.get_res import *

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot the corner plot of Df with theta_E of the main lens")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    
    ####################
    present_program(sys.argv[0])
    ####################
    
    labels = [r"$\theta^1_E$","$\Delta\phi_{Fermat} AB$", "$\Delta\phi_{Fermat} AC$","$\Delta\phi_{Fermat} BC$"]

    for sets in settings:
        savefig_path = get_savefigpath(sets)
        f = get_mcmc_fermat(sets).T
        df = np.array([f[1]-f[0],f[2]-f[0],f[2]-f[1]]) #AB,AC,BC
        smpl  = get_mcmc_smpl(sets)
        param = get_mcmc_prm(sets)
        for i,prm in enumerate(param):
            if "theta_E_lens0" in prm:
                theta_E = smpl.T[i]
        crn = [theta_E,*df.tolist()]
        title = sets.replace("_"," ").replace(".py","").replace("settings","")
        plot = corner.corner(np.transpose(crn), labels=labels, show_titles=True)
        plot.suptitle(title, fontsize=16)
        plot.savefig(savefig_path+"/Corner_Df_thetaE.pdf")
    success(sys.argv[0])

