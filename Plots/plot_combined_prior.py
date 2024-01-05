#!/usr/bin/env python
# coding: utf-8

import copy
import json
import os,sys 
import argparse
import numpy as np

from Utils.tools import *
from Utils.get_res import load_whatever
from Prior.Prior import Df_prior, Df_prior_ABC
from Plots.plotting_tools import plot_probability3D
from Utils.statistical_tools import get_bins_volume

def get_minmax(samples,dim):
    # We get the min and maximum for each dimension relatively to all the datasets togheter
    return np.min([min(par_i) for par_i in np.array(samples,dtype=object)[:,dim]]),\
            np.max([max(par_i) for par_i in np.array(samples,dtype=object)[:,dim]])

def KDE_prior(Priors,d_bins):
    # Priors shape must be (filter_i, point_i, dimensions), eg: (1,1000,3)
    # (same as standard smpl)
    # return Priors evaluated with KDE at center of the bins
    
    from sklearn.neighbors import KernelDensity
    from sklearn.model_selection import GridSearchCV
    
    # find the best bandwidth for each filter
    prior_bndws = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth': prior_bndws})
    # initial shape of samples: filter,dimension,step

    #bandwidths = [ grid.fit(np.transpose(samples[filt_i])) for filt_i in range(len(filters)) ]   
    bandwidths = [grid.fit(Prior_i) for Prior_i in Priors]   
    # input of grid.fit must have shape: ( n*steps, n*dimensions)

    # define the KDE
    kdes = [KernelDensity(**bndw.best_params_, kernel='gaussian') for bndw in bandwidths]
    # Fit the KDEs to the samples
    KDEs = [kde.fit(Priors[i]) for i,kde in enumerate(kdes)]

    # we then have to sample the space
    positions = KDE_position_from_bin(d_bins) #position.shape = ((n_bins-1)**Dim,Dim) 
    nbins  = [len(b)-1 for b in d_bins]
    Priors = [kde_i.score_samples(positions).reshape(nbins[0],nbins[1],nbins[2]) for kde_i in KDEs]
    
    # NOTE: by construction the KDE is not normalised
    # we could "normalise" it pointwise, but it is actually not necessary, 
    # as it is only shifted by a factor, and the Combined_PDF has to be normalised
    # anyway after the multiplication
    # but this is blocking the test so
    Priors = Priors/np.sum(Priors) # normalised pointwise
    return Priors, positions

def KDE_position_from_bin(d_bins):
    # obtain center of 3D bins and return as
    # positions for KDE -> shape: ((nbins-1)**Dim,Dim)
    if len(d_bins)!=3:
        raise ValueError("Implemented explicitely for 3 dimensions")

    # get the center of the bins instead of the edges
    centers = []
    for bi in d_bins:
        cnt_di = [] # for each dim_i indep
        for i in range(len(bi)-1):
            cnt_di.append(.5*(bi[i]+bi[i+1]))
        centers.append(cnt_di)
    # From here on considering 3Dims
    #
    # Creating a grid of points from the matrix of centers
    Dim1, Dim2, Dim3 = np.meshgrid(*centers)
    # reorganising the shape of it
    positions = np.vstack([Dim1.ravel(), Dim2.ravel(), Dim3.ravel()]).T    
    return positions
    
def check_combined(Combined_PDF):
    #Check superposition
    if np.sum(Combined_PDF)==0.0:
        raise RuntimeError("No superposition between posteriors. No idea yet how to solve it")


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot the multiplied posterior distribution of the Fermat potential difference from the given filters")

    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                        help="Directory name where to save the multiplied posteriors")
    parser.add_argument("-nb","--number_bins",type=int, dest="nbins", default=100,
                        help="Number of bins per dimension (be careful with it! too many bins can be catastrophic)")
    parser.add_argument("-mp","--mcmc_prior",type=int, dest="mcmc_prior", default=1000,
                        help="Number of steps for the MCMC sampling of the Priors")
    parser.add_argument("-KDE", action="store_true", dest="KDE", default=False,
                        help="Use KDE (Kernel Density Estimator) instead of histograms (WARNING:Very slow for high number of points and/or bins)")
    parser.add_argument("-mcmc","--MCMC", action="store_true", dest="mcmc", default=False,
                        help="Also do the MCMC integration of the posterior")

    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

    args = parser.parse_args()
    dir_name = args.dir_name
    nbins = int(args.nbins)
    mcmc  = bool(args.mcmc)
    mcmc_prior = int(args.mcmc_prior)
    cut_mcmc = 0
    setting_names =  args.SETTING_FILES  
    ####################
    present_program(sys.argv[0])
    ####################

    setting_name = []
    setting_position = []
    for sn in setting_names:
        setting_position.append(find_setting_path(sn))
        setting_name.append(sn.replace(".py",""))

    savemcmc_path=[]
    for i in setting_name:
        path_i="./backup_results/"+i.replace("settings_","mcmc_")+"/"
        savemcmc_path.append(path_i)


    # In[ ]:



    fermat_mcmc=[]
    for i in range(len(setting_name)):
        fermat_file = savemcmc_path[i]+setting_name[i].replace("settings","mcmc_ordered_fermat")+".json"
        if not os.path.isfile(fermat_file):
            create_fermat = "python fermat_pot_analysis.py "+setting_name[i]
            os.system(create_fermat)
        fermat_mcmc.append(fermat_file)

    samples = []
    for i in fermat_mcmc:
        with open(i, 'r') as mcmc_file_i:
            mcmc_iT = np.array(json.load(mcmc_file_i))
        cut_mcmc_scaled = int(len(mcmc_iT)*cut_mcmc/1000)
        mcmc_i   = np.transpose(mcmc_iT[cut_mcmc_scaled:]) 
        mcmc_Dfi = mcmc_i[1:]-mcmc_i[0]
        #############################################
        #############################################    
        mcmc_c    = np.array(copy.deepcopy(mcmc_Dfi))
        mcmc_BC   = mcmc_c[1] - mcmc_c[0]  # BC = C - B = (C-A)-(B-A) = AC - AB
        mcmc_c[2] = mcmc_BC
        mcmc_Dfi  = mcmc_c 
        samples.append(mcmc_Dfi.tolist())
    Df_boundaries   = [get_minmax(samples,i) for i in range(3)]

    Priors = [Df_prior_ABC(sett,npoints=mcmc_prior,Df_boundaries=Df_boundaries) for sett in setting_names] #those are MCMC priors

    N = len(samples)
    param_dim = len(samples[0]) # the parametric dimension must be the same for all datasets 
    #First find the same bins for each dataset in all parametric dimensions
    DfABmin,DfABmax = Df_boundaries[0]
    DfACmin,DfACmax = Df_boundaries[1]
    DfADmin,DfADmax = Df_boundaries[2]
    maxs,mins       = [DfABmax,DfACmax,DfADmax],[DfABmin,DfACmin,DfADmin]
    d_bins =[]
    for par_dim_i in range(param_dim): # for each parametric dim
        # Find for each dataset the corresponding max and min
        max_dimi = maxs[par_dim_i]
        min_dimi = mins[par_dim_i]
        diff_dimi = max_dimi-min_dimi
        len_bin_dimi = diff_dimi/nbins
        d_bins_dimi = [min_dimi]
        for nbin_i in range(nbins):
            d_bins_dimi.append(d_bins_dimi[-1]+len_bin_dimi)
        d_bins.append(d_bins_dimi)
        
        
    Priors,Positions_KDE_Prior = KDE_prior(Priors,d_bins)


    Combined_Prior = Priors[0]
    for i in range(1,N):
        Combined_Prior *= Priors[i]
    # Check that the Combined Prior is not 0
    check_combined(Combined_Prior)


    save_dir= "./backup_results/PDF_multiplication_ABC"
    mkdir(save_dir)


    filters=[]
    for i in setting_name:
        filters.append(i.replace("settings_",""))
        
    if dir_name==".":
        lf = ""
        for f in filters:
            lf+=f+"_"
        lf = lf[:-1]
        save_dir=save_dir/str(lf)
    else:
        save_dir=save_dir/dir_name
    mkdir(save_dir)

    Combined_bins = load_whatever(save_dir/str("Combined_PDF_bins.pkl") )
    
    ########
    print("WARNING: This is a test in ",sys.argv[0]," to see if the Combined_prior and the original prior are similar enough")


    # not sure this is legit:
    Combined_Prior/= np.sum(Combined_Prior*get_bins_volume(Combined_bins))
    Priors_i     = Priors[0]/np.sum(Priors[0]*get_bins_volume(Combined_bins))
    plot = plot_probability3D(Combined_Prior,Combined_bins,labels=["AB","AC","AD"],udm="\"",alpha=0.3)
    plot.savefig(str(save_dir)+"/Compare_prior_comb.png")
    plot = plot_probability3D(Priors_i,Combined_bins,labels=["AB","AC","AD"],colors=["y","k","b"],udm="\"",alpha=0.3)
    plot.savefig(str(save_dir)+"/Compare_prior_1.png")
    print("WARNING: checkout "+str(save_dir)+"/tmp_Compare_prior.png\n")

    success(sys.argv[0])
