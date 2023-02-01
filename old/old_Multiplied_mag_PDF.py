#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Copy from Multiply_Df_PDF
# now considering multiplication of mag posteriors


# In[ ]:


import numpy as np
import json
import os
import matplotlib.pyplot as plt
import pickle
from corner import corner
import sys
import pathlib as pth
import argparse
from lenstronomy.Sampling.Pool.multiprocessing import MultiPool
from lenstronomy.Util.sampling_util import sample_ball
#from mag_remastered import mag
from emcee import EnsembleSampler
import copy

#my libs
from Prior import mag_prior_ABC
from mag_remastered import labels
from Multiply_PDF import Multiply_PDF,get_minmax
from tools import *
from statistical_tools import get_bins_volume


# In[ ]:


parser = argparse.ArgumentParser(description="Plot the multiplied posterior distribution of magnification ratio from the given filters")
parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                    help="Cut the first <c> steps of the mcmc to ignore them")
parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the multiplied posteriors")
parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                    help="Number of bins per dimension (be careful with it! too many bins can be catastrophic)")
parser.add_argument("-ms","--mcmc_steps",type=int, dest="mcmc_steps", default=1000,
                    help="Number of steps for the MCMC sampling and plot")
parser.add_argument("-mp","--mcmc_prior",type=int, dest="mcmc_prior", default=1000,
                    help="Number of steps for the MCMC sampling of the Priors")
parser.add_argument("-KDE", action="store_true", dest="KDE", default=False,
                    help="Use KDE (Kernel Density Estimator) instead of histograms (WARNING:Very slow for high number of points and/or bins)")
parser.add_argument("-mcmc","--MCMC", action="store_true", dest="mcmc", default=False,
                    help="Also do the MCMC integration of the posterior")
parser.add_argument("-b","--boundaries", action="store_true", dest="boundaries", default=False,
                    help="Consider more precise but long prior: computer mag boundaries and sample only there")

parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

args = parser.parse_args()
cut_mcmc = int(args.cut_mcmc)
dir_name = args.dir_name
KDE = bool(args.KDE)
nbins = int(args.nbins)
mcmc  = bool(args.mcmc)
mcmc_steps  = int(args.mcmc_steps)
mcmc_prior  = int(args.mcmc_prior)
boundaries  = bool(args.boundaries)

setting_names =  args.SETTING_FILES  
####################
present_program(sys.argv[0])
####################

setting_name = []
setting_position = []
for sn in setting_names:
    setting_position.append(find_setting_path(sn))
    setting_name.append(sn.replace(".py",""))

backup_path="backup_results"

filters = [strip_setting_name(st) for st in in setting_name]

save_dir = create_dir_name(setting_names,save_dir="PDF_multiplication_ABC",dir_name=dir_name,backup_path=backup_path)

savemcmc_path = [get_savemcmcpath(st) for st in  setting_name]


for st_i in range(len(setting_position)):
    os.system("cp "+str(setting_position[st_i])+"/"+str(setting_names[st_i])+\
    " "+str(save_dir)+"/.")


# for mag ratio
samples = []
for i in range(len(setting_name)):
    mcmc_iT  = get_mcmc_mag(setting_name[i]) 
    cut_mcmc_scaled = int(len(mcmc_iT)*cut_mcmc/1000)
    mcmc_mag_i   = np.transpose(mcmc_iT[cut_mcmc_scaled:]) # shape: (dimensions,steps)
    #############################################
    #############################################    
    #print("WARNING: Given the low S/N of image D, I will discard here and instead consider Mag BC")    
    mcmc_c     = np.array(copy.deepcopy(mcmc_mag_i))
    mcmc_BC    = mcmc_c[1]/mcmc_c[0]  # BC = C / B = (C/A)/(B/A) = AC/AB
    mcmc_c[2]  = mcmc_BC
    mcmc_mag_i = mcmc_c 
    samples.append(mcmc_mag_i.tolist())
print("WARNING: Given the low S/N of image D, I will discard here and instead consider Mag BC")    
param_names = [*labels[:3],"$\mu_C$/$\mu_B$"]


# In[ ]:


# If using it on my machine with dark style jupyter
if my_machine():
    plt.style.use(['classic'])


# In[ ]:


# First we must compute the prior for the mag ratio P(mu)

mag_boundaries = None
if boundaries:
    mag_boundaries = [get_minmax(samples,i) for i in range(3)]                   
Priors_mag = [mag_prior_ABC(sett,save_mcmc=True,npoints=mcmc_prior,mag_boundaries=mag_boundaries) for sett in setting_names] #those are MCMC priors
pt = corner(np.array(Priors_mag[0]),labels=param_names,show_titles=True)
pt.savefig(str(save_dir)+"/mcmc_prior_mag_"+setting_names[0].replace(".py","")+".pdf")

# NOTE: not the same shape as the samples!
# samples.shape = (n filter, n dim, n point)*,    Priors_mag.shape = (n filter, n points, n dims)
# *: n point is usually not outputed by the np.shape bc samples have different number of points

#These priors are then multiplied and combined into a single prior using the same binning as
# the Combined_PDF so that it is the comparable and they are then combined according to 
# the theory (see Notes.ipynb at the 16th of Feb)

if KDE:
    from Multiply_PDF import Multiply_PDF_KDE
    npoints = nbins  # To TEST
    Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Priors=Priors_Df,savedir=save_dir)
else:
    from Multiply_PDF import Multiply_PDF_HIST
    Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Priors=Priors_Df,savedir=save_dir)

# The Combined_PDF must be re-normalised
if KDE:
    #Combined_PDF = Normalise_KDE(Combined_PDF,Positions_KDE)
    Combined_PDF /= np.sum(Combined_PDF)
    # still unclear how to do it in KDE
else:
    Combined_PDF/=np.sum(Combined_PDF*get_bins_volume(Combined_bins))


with open(str(save_dir)+"/Combined_mag_PDF"+["_KDE" if KDE else ""][0]+".pkl","wb") as f:
    pickle.dump(Combined_PDF,f)
    
if KDE:
    with open(str(save_dir)+"/Combined_mag_PDF_KDE_positions.pkl","wb") as f:
        pickle.dump(positions,f)
else:
    with open(str(save_dir)+"/Combined_mag_PDF_bins.pkl","wb") as f:
        pickle.dump(Combined_bins,f)


# In[ ]:


# We need to sample it for the plot 
if not KDE and mcmc:
    from statistical_tools import *
    from multiprocessing import cpu_count
    mcmc_init_pos,mcmc_sigma = estimate_for_mcmc(Combined_PDF,Combined_bins)
    #mcmc_sampling = sampler(mcmc_init_pos,prob=Combined_PDF,bins=Combined_bins,mcmc_simga=mcmc_sigma,mcmc_steps=int(1e6))
    #mcmc_chain = mcmc_sampling[0]
    def logP(pos):
        prob_at_pos = get_prob_at_pos(pos,Combined_PDF,Combined_bins)
        if prob_at_pos==0:
            return -np.inf
        return np.log(prob_at_pos)

    nprocesses = cpu_count()-1
    nwalkers = 42
    init_sample = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
    pool = MultiPool(processes=nprocesses) 
    emcee_mcmc = EnsembleSampler(nwalkers,len(np.shape(Combined_PDF)), logP,pool=pool)
    initial_state = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
    emcee_mcmc.run_mcmc(initial_state= initial_state,nsteps=mcmc_steps)
    with open(str(save_dir)+"/mcmc_mag_chain.json","w") as f:
        json.dump(np.array(mcmc_chain).tolist(),f)
    plot = corner(mcmc_chain,bins=nbins,labels=[ p+" [\"]" for p in param_names], show_titles=True)
    plot.savefig(str(save_dir)+"/MCMC_multiplication_mag.png")

if not KDE:
    from plotting_tools import plot_probability3D
    plot = plot_probability3D(Combined_PDF,Combined_bins,labels=param_names,udm="\"")
    plot.savefig(str(save_dir)+"/CombinedProbability_mag.pdf")
else:
    from plotting_tools import plot_probability3D_KDE
    plot = plot_probability3D_KDE(Combined_PDF,Positions_KDE,labels=param_names,udm="\"")
    plot.savefig(str(save_dir)+"/CombinedProbability_KDE_mag.pdf")


# In[ ]:


print("Result directory:", str(save_dir))
success(sys.argv[0])

