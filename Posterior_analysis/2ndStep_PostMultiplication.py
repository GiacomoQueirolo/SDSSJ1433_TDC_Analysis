#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Idea: use the Multiplied_Df_PDF program to define in which boundary of the parametric space the posterior are overlapping, 
# then compute a kde only within that region for each filter posterior with only the subset of point contained in there
# and then multiply them togheter

# the Prior correction shows it works in that case:


# In[2]:


import copy
import os,sys
import argparse
import numpy as np
import json,pickle
from os import walk
import pathlib as pth
from corner import corner
import matplotlib.pyplot as plt

#my libs
from Utils.tools import *
from Utils.get_res import *
from Posterior_analysis.Multiply_PDF import get_minmax
from Prior import Df_prior, Df_prior_ABC
from Utils.statistical_tools import get_bins_volume
from Utils.combinedsetting_class import combined_setting


# In[ ]:


parser = argparse.ArgumentParser(description="Plot the multiplied posterior distribution of the Fermat potential difference from the given filters - 2nd step",
                            usage='pokus --help',)
parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the multiplied posteriors")
"""parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=40,
                    help="Number of bins per dimension (be careful with it! too many bins can be catastrophic)")
parser.add_argument("-ms","--mcmc_steps",type=int, dest="mcmc_steps", default=1000,
                    help="Number of steps for the MCMC sampling and plot")
parser.add_argument("-mp","--mcmc_prior",type=int, dest="mcmc_prior", default=1000,
                    help="Number of steps for the MCMC sampling of the Priors")
parser.add_argument("-KDE", action="store_true", dest="KDE", default=False,
                    help="Use KDE (Kernel Density Estimator) instead of histograms (WARNING:Very slow for high number of points and/or bins)")
parser.add_argument("-mcmc","--MCMC", action="store_true", dest="mcmc", default=False,
                    help="Also do the MCMC integration of the posterior")
#parser.add_argument("-b","--boundaries", action="store_true", dest="boundaries", default=False,
#                    help="Consider more precise but long prior: computer Df boundaries and sample only there")
parser.add_argument("-NOP","--not_old_Prior", action="store_false", dest="old_prior", default=True,
                    help="If present it will compute the Prior again (might take a while), else look for previously computed ones")
"""
parser.add_argument('POST_DIR',nargs="+",default=[],help="Path to the directory of the previously computed combined posterior using Multiplied_Df_PDF.py with histograms")

args = parser.parse_args()
dir_name = args.dir_name
"""cut_mcmc = int(args.cut_mcmc)
KDE = bool(args.KDE)
nbins = int(args.nbins)
mcmc  = bool(args.mcmc)
mcmc_steps  = int(args.mcmc_steps)
mcmc_prior  = int(args.mcmc_prior)
old_prior   = bool(args.old_prior)
setting_names =  args.SETTING_FILES  
"""
data_path = args.POST_DIR
backup_path = "backup_results/"
if backup_path not in data_path:
    data_path = backup_path+data_path
list_data = next(walk(data_path), (None, None, []))[2]
settings = []
for ld in list_data :
    if "setting" in ld:
        settings.append(ls)
combined_prob = data_path+"/"

