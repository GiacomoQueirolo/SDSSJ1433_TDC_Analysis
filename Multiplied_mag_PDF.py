#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Copy from Multiply_Df_PDF 22nd Dec '22

import copy
import os,sys
import argparse
import numpy as np
import json,pickle
import pathlib as pth
from corner import corner
import matplotlib.pyplot as plt
from lenstronomy.Util.sampling_util import sample_ball
from lenstronomy.Sampling.Pool.multiprocessing import MultiPool

#my libs
from tools import *
from get_res import *
from Prior import mag_prior_ABC
from mag_remastered import labels
from Multiply_PDF import get_minmax
from statistical_tools import get_bins_volume


# In[ ]:


parser = argparse.ArgumentParser(description="Plot the multiplied posterior distribution of the magnification ratio from the given filters",
                            usage='pokus --help',)
parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                    help="Cut the first <c> steps of the mcmc to ignore them")
parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the multiplied posteriors")
parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=40,
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

parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

args = parser.parse_args()
cut_mcmc = int(args.cut_mcmc)
dir_name = args.dir_name
KDE = bool(args.KDE)
nbins = int(args.nbins)
mcmc  = bool(args.mcmc)
mcmc_steps  = int(args.mcmc_steps)
mcmc_prior  = int(args.mcmc_prior)
old_prior   = bool(args.old_prior)
setting_names =  args.SETTING_FILES  
####################
present_program(sys.argv[0])
####################

backup_path ="backup_results"
main_savedir="PDF_multiplication_ABC"

filters       = [get_filter(st) for st in setting_names]
savemcmc_path = [get_savemcmcpath(st) for st in  setting_names]
save_dir  = create_dir_name(setting_names,save_dir=main_savedir,dir_name="Mag",backup_path=backup_path,copy_settings=False)
save_dir  = create_dir_name(setting_names,save_dir=save_dir.replace(backup_path,""),dir_name=dir_name,backup_path=backup_path,copy_settings=True)

save_log_command(save_dir)
# for mag ratio
samples = []
for i in range(len(setting_names)):
    mcmc_iT  = get_mcmc_mag(setting_names[i]) 
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
labels = [*labels[:3],"$\mu_C$/$\mu_B$"] 
##################################################################################


# In[ ]:


# If using it on my machine with dark style jupyter
if my_machine():
    plt.style.use(['classic'])


# In[ ]:


# Then we make the 3D histogram #-> should i consider more complex idea? like kernels?
# actually Kernels sounds better -> see KDE_Df_PDF and old/KDE_Df_PDF_example
# the question then is, how to multiply them? 
# fit them at the bin center and multiply there -> another idea would be to keep the various KDE and compute them
# only when needed at the needed position/ at no matter which binning size

# we could do that But keep it going with the multiplication at the bin position for the plots


# In[ ]:


#Try to find the best one 
# -<> sklearn I know and it works, the other seems easier and more reliable
# following https://jakevdp.github.io/PythonDataScienceHandbook/05.13-kernel-density-estimation.html
# i go for sklearn: "While there are several versions of kernel density estimation
# implemented in Python (notably in the SciPy and StatsModels packages),
# I prefer to use Scikit-Learn's version because of its efficiency and flexibility.  "
# see Kernel-Density-Estimation_new.ipynb for the full notebook
# SEE OLD NOTES IN old/KDE_exercise_2Feb.py


# In[ ]:


# First we must compute the prior for the fermat potential P(Delta Phi)

#Priors_Df = [Df_prior(setting_names[0],npoints=mcmc_prior)] #those are MCMC priors
#name_prior_HRes = "Prior_"+"_".join(filters)+"_highres.json"
#path_prior_HRes = backup_path+"/"+main_savedir+"/"+name_prior_HRes
name_prior = "kw_Prior_mag_"+"_".join(filters)+".json"
path_prior = backup_path+"/"+main_savedir+"/"+name_prior
#name_prior_LRes = "Prior_"+"_".join(filters)+"_lowres.json"
#path_prior_LRes = backup_path+"/"+main_savedir+"/"+name_prior_LRes
#no_old_prior_hr = False
no_old_prior = False
if old_prior:
    try:
        #Prior_Df_HRes = load_whatever(path_prior_HRes)
        kw_Prior_mag = load_whatever(path_prior)
        Prior_mag    = kw_Prior_mag["Prior"]
        NtotPrior    = kw_Prior_mag["Ntot"]
    except:
        #no_old_prior_hr = True
        no_old_prior = True
    """try:
        Prior_Df_LRes = load_whatever(path_prior_LRes)
    except:
        no_old_prior_lr = True
    """
#if no_old_prior_hr or no_old_prior_lr or not old_prior:
if no_old_prior or not old_prior:
    mag_boundaries = [get_minmax(samples,i) for i in range(3)]
    
    #if no_old_prior_hr or not old_prior:
    if no_old_prior or not old_prior:
        """
        Df_boundaries = None
        if boundaries:
            Df_boundaries = [get_minmax(samples,i) for i in range(3)]                   
        Priors_Df = [Df_prior_ABC(sett,save_mcmc=True,npoints=min(mcmc_prior,500000),Df_boundaries=Df_boundaries) for sett in setting_names] #those are MCMC priors
        """
        #Prior_Df = Df_prior_ABC(setting_names[0],save_mcmc=True,npoints=min(mcmc_prior,500000),Df_boundaries=Df_boundaries) 
        # We actually assume them to be ~identical priors 
        # we consider a high res. prior on the bounded region and a 
        # lower reson a larger region
        #Prior_Df_HRes,NtotPrior = Df_prior_ABC(setting_names[0],save_mcmc=True,npoints=min(mcmc_prior,50000),Df_boundaries=Df_boundaries) 
        Prior_mag,NtotPrior = mag_prior_ABC(setting_names[0],save_mcmc=True,npoints=min(mcmc_prior,50000),mag_boundaries=mag_boundaries) 
        #Prior_Df_LRes = Df_prior_ABC(setting_names[0],save_mcmc=True,npoints=min(mcmc_prior,50000),Df_boundaries=None,output_name="mcmc_prior_Df_lowres")
        #pt = corner(np.array(Priors_Df[0]),labels=param_names,show_titles=True)
        #pt.savefig(str(save_dir)+"/mcmc_prior_Df_"+setting_names[0].replace(".py","")+".pdf")
        #pt = corner(np.array(Prior_Df_LRes),labels=param_names,show_titles=True)
        #pt.savefig(str(save_dir)+"/mcmc_prior_Df_lowres"+setting_names[0].replace(".py","")+".pdf")
        #pt = corner(np.array(Prior_Df_HRes),labels=param_names,show_titles=True)
        pt = corner(np.array(Prior_mag),labels=labels,show_titles=True)
        #pt.savefig(str(save_dir)+"/mcmc_prior_Df_highres"+setting_names[0].replace(".py","")+".pdf")
        pt.savefig(str(save_dir)+"/mcmc_prior_mag"+setting_names[0].replace(".py","")+".pdf")
        # NOTE: not the same shape as the samples!
        # samples.shape = (n filter, n dim, n point)*,    Priors_Df.shape = (n filter, n points, n dims)
        # *: n point is usually not outputed by the np.shape bc samples have different number of points

            #These priors are then multiplied and combined into a single prior using the same binning as
            # the Combined_PDF so that it is the comparable and they are then combined according to 
            # the theory (see Notes.ipynb at the 16th of Feb)
            # -> actually is only 1 prior, bc we assume them to be compatible
        #save_json(Prior_Df_HRes,path_prior_HRes)
        kw_Prior_mag = {"Prior":np.array(Prior_mag).tolist(),"Ntot":NtotPrior}
        save_json(kw_Prior_mag,path_prior)
        #save_json(Prior_Df_LRes,path_prior_LRes)
"""    if no_old_prior_lr or not old_prior:
        Prior_Df_LRes = Df_prior_ABC(setting_names[0],save_mcmc=True,npoints=min(mcmc_prior,50000),Df_boundaries=None,output_name="mcmc_prior_Df_lowres")
        pt = corner(np.array(Prior_Df_LRes),labels=param_names,show_titles=True)
        pt.savefig(str(save_dir)+"/mcmc_prior_Df_lowres"+setting_names[0].replace(".py","")+".pdf")
        save_json(Prior_Df_LRes,path_prior_LRes)
"""
if KDE:
    from Multiply_PDF import Multiply_PDF_KDE
    npoints = nbins  # To TEST
    #Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Priors=Priors_Df,savedir=save_dir)
    #Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Prior_LR=Prior_Df_LRes,Prior_HR=Prior_Df_HRes,savedir=save_dir)
    Combined_PDF,Positions_KDE = Multiply_PDF_KDE(samples,npoints,Prior=Prior_mag,NtotPrior=NtotPrior,savedir=save_dir)
else:
    from Multiply_PDF import Multiply_PDF_HIST
    #Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Priors=Priors_Df,savedir=save_dir)
    #Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Prior_LR=Prior_Df_LRes,Prior_HR=Prior_Df_HRes,savedir=save_dir)
    Combined_PDF,Combined_bins = Multiply_PDF_HIST(samples,nbins,Prior=Prior_mag,NtotPrior=NtotPrior,savedir=save_dir)



# The Combined_PDF must be re-normalised
if KDE:
    #Combined_PDF = Normalise_KDE(Combined_PDF,Positions_KDE)
    Combined_PDF /= np.sum(Combined_PDF)
    # still unclear how to do it in KDE
else:
    Combined_PDF /= np.sum(Combined_PDF*get_bins_volume(Combined_bins))

with open(str(save_dir)+"/Combined_mag_PDF"+["_KDE" if KDE else ""][0]+".pkl","wb") as f:
    pickle.dump(Combined_PDF,f)
    
if KDE:
    with open(str(save_dir)+"/Combined_mag_PDF_KDE_positions.pkl","wb") as f:
        pickle.dump(positions,f)
else:
    with open(str(save_dir)+"/Combined_mag_PDF_bins.pkl","wb") as f:
        pickle.dump(Combined_bins,f)



from Multiplied_Df_PDF import combined_setting

comment="Product of posteriors done by "+str(sys.argv[0])
z_lens = [s.z_lens for s in get_setting_module(setting_names,1) ]
z_source = [s.z_source for s in get_setting_module(setting_names,1) ]
if all(zl==z_lens[0] for zl in z_lens) and all(zs==z_source[0] for zs in z_source):
    z_lens = z_lens[0]
    z_source = z_source[0]
else:
    raise RuntimeError("Not all settings have same redshift for lens and source. Something is off")
CombSett=combined_setting(comment,z_source,z_lens,filters,setting_names)
pickle.dump(CombSett,open(str(save_dir)+"/combined_setting.pkl","wb"))
#####################################################################################

# We need to sample it for the plot 
if not KDE and mcmc:
    from statistical_tools import *
    from emcee import EnsembleSampler
    from multiprocessing import cpu_count
    mcmc_init_pos,mcmc_sigma = estimate_for_mcmc(Combined_PDF,Combined_bins)
    #mcmc_sampling = sampler(mcmc_init_pos,prob=Combined_PDF,bins=Combined_bins,mcmc_simga=mcmc_sigma,mcmc_steps=int(1e6))
    #mcmc_chain = mcmc_sampling[0]
    def logP(pos):
        prob_at_pos = get_prob_at_pos(pos,Combined_PDF,Combined_bins)
        if prob_at_pos==0:
            return -np.inf
        return np.log(prob_at_pos)

    nprocesses    = cpu_count()-1
    nwalkers      = 42
    init_sample   = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
    pool          = MultiPool(processes=nprocesses) 
    emcee_mcmc    = EnsembleSampler(nwalkers,len(np.shape(Combined_PDF)), logP,pool=pool)
    initial_state = sample_ball(mcmc_init_pos,mcmc_sigma,nwalkers)
    emcee_mcmc.run_mcmc(initial_state= initial_state,nsteps=mcmc_steps)
    mcmc_chain = emcee_mcmc.get_chain(discard=0, thin=1, flat=True)
    with open(str(save_dir)+"/mcmc_mag_chain.json","w") as f:
        json.dump(np.array(mcmc_chain).tolist(),f)
    plot = corner(mcmc_chain,bins=nbins,labels=[ p+" []" for p in labels], show_titles=True)
    plot.savefig(str(save_dir)+"/MCMC_multiplication_mag.png")

if not KDE:
    from plotting_tools import plot_probability3D
    plot = plot_probability3D(Combined_PDF,Combined_bins,labels=labels,udm="")
    plot.savefig(str(save_dir)+"/CombinedProbability_mag.pdf")
else:
    from plotting_tools import plot_probability3D_KDE
    plot = plot_probability3D_KDE(Combined_PDF,Positions_KDE,labels=labels,udm="")
    plot.savefig(str(save_dir)+"/CombinedProbability_KDE_mag.pdf")



# In[ ]:


print("Result directory:", str(save_dir))
success(sys.argv[0])

