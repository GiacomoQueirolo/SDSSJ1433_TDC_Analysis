#!/usr/bin/env python
# coding: utf-8



import pickle
import numpy as np
import scipy.stats as st
from multiprocessing import Pool
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV


from Utils.statistical_tools import get_bins_volume
from Plots.plotting_tools import plot_probability3D,plot_probability3D_KDE

labels_BC,labels_AD = ["AB","AC","BC"],["AB","AC","AD"]


def get_minmax(samples,dim=None):
    # We get the min and maximum for each dimension relatively to all the datasets togheter
    if dim is None:
        return np.min([np.min(s,axis=1) for s in samples],axis=0), np.max([np.max(s,axis=1) for s in samples],axis=0)
    else:
        return np.min([np.min(s,axis=1) for s in samples],axis=0)[dim], np.max([np.max(s,axis=1) for s in samples],axis=0)[dim]


def Multiply_PDF_HIST_fitPrior(samples,nbins,Prior,savedir=".",labels=labels_AD,verbose=True):
    # Using histograms to obtain the binned density for each datasets
    N = len(samples)
    Dim = len(samples[0])
    #First find the same bins for each dataset in all parametric dimensions
    mins,maxs  = get_minmax(samples)#[DfABmax,DfACmax,DfBCmax],[DfABmin,DfACmin,DfBCmin]
    len_bins = (maxs-mins)/nbins
    d_bins = np.transpose([mins +i*len_bins for i in range(nbins+1)])
    
    # Then we make the histogram 
    PDFs,Combined_bins = hists_given_bins(samples,d_bins)
    

    Prior,_ = KDE_fit_prior(Prior=Prior,d_bins=d_bins,verbose=verbose)
            
    Combined_PDF = PDFs[0]
    for i in range(1,N):
        Combined_PDF *= PDFs[i]
    
    # Check that the Combined PDF is not 0
    check_combined(Combined_PDF,PDFs,Combined_bins,HIST="True")

    if Prior is not None:
        check_combined(Prior,KDE=True)
        Norm_Prior = Prior/np.sum(Prior*get_bins_volume(d_bins))
        if Dim==3:
            plot = plot_probability3D(Norm_Prior,Combined_bins,labels=labels,udm="\"",alpha=0.3)
            plot.savefig(str(savedir)+"/Normalised_prior.png")
        for i in range(N-1):
            Combined_PDF/=Norm_Prior                
    return Combined_PDF, Combined_bins


def check_combined(Combined_PDF,PDFs=None,PDF_bins=None,HIST=False,KDE=False):
    #Check superposition
    if HIST:
        if np.sum(Combined_PDF)==0.0:
            for i in range(len(PDFs)):
                print(i," Max :",PDF_bins[np.where(PDFs[i]==np.max(PDFs[i]))[0]])
            raise RuntimeError("No superposition between posteriors. \
            Try adjusting the bin nummer or correcting the modelling procedure for one or more dataset")
    elif KDE:
        if np.sum(Combined_PDF)==0.0:
            raise RuntimeError("No superposition between posteriors. No idea yet how to solve it")
    else:
        raise RuntimeError("Give either KDE or HIST")

def hists_given_bins(samples,bins):
    PDFs = []
    PDFs_bins = []
    for i in range(len(samples)):
        smpl = np.array(samples[i]).tolist()
        density,bins = np.histogramdd(smpl,bins=bins,density=True) 
        # density is norm given the area/volume of the bin
        PDFs.append(density)
        PDFs_bins.append(bins)
        
    #Sanity check
    same_binning = True
    for i in range(1,len(PDFs_bins)):
        if np.all(PDFs_bins[i]!=PDFs_bins[0]):
            same_binning = False
            raise ValueError("We have a problem, we should have the same binning for each filter")
    if same_binning:
        PDF_bins = PDFs_bins[0]
        
    return PDFs,PDF_bins


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


def KDE_fit_prior(Prior,d_bins,verbose=True):
    # Prior is given as class 
    # not uniform! as it's actually the prior in the fermat space
    # BUT sampled all over it, not only the region of interest
     # return Priors evaluated with KDE at center of the bins
    
    # find the best bandwidth for each filter
    prior_bndws = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth': prior_bndws})
    # initial shape of samples: filter,dimension,step
    prior_Df_sample = Prior.get_Df_sample()
    bandwidth = grid.fit(prior_Df_sample)
    # define the KDE
    if verbose:
        print("Defining KDE")
    kde = KernelDensity(**bandwidth.best_params_, kernel='gaussian') 
    
    # Fit the KDE to the samples
    if verbose:
        print("Fitting KDE to the prior sample")
    KDE = kde.fit(prior_Df_sample)
    
    # we then have to sample the space
    positions = KDE_position_from_bin(d_bins) #position.shape = ((n_bins-1)**Dim,Dim) 
    shape     = [len(b)-1 for b in d_bins]
    if verbose:
        print("Getting the value at the desired positions")
    Prior     = KDE.score_samples(positions).reshape(*shape)
    return Prior, positions


############################################################
####################### OLD and WIP: ####################### 
############################################################
    
def Multiply_PDF_HIST(samples,nbins,Prior=None,NtotPrior=None,savedir="."):
    # Using histograms to obtain the 3D binned density for each datasets
    N = len(samples)
    Dim = len(samples[0])
    #First find the same bins for each dataset in all parametric dimensions
    mins,maxs  = get_minmax(samples)#[DfABmax,DfACmax,DfBCmax],[DfABmin,DfACmin,DfBCmin]
    len_bins = (maxs-mins)/nbins
    d_bins = np.transpose([mins +i*len_bins for i in range(nbins+1)])
    
    # Then we make the 3D histogram    
    PDFs,Combined_bins = hists_given_bins(samples,d_bins)
    
    if Prior is not None: # and Prior_LR is not None
        #Prior,Positions_KDE_Prior = KDE_prior_corr(Prior_HR=Prior_HR,Prior_LR=Prior_LR,d_bins=d_bins)
        Prior,Positions_KDE_Prior = KDE_prior_corr(Prior=Prior,NtotPrior=NtotPrior,d_bins=d_bins)
        #Priors,Prior_bins = hists_given_bins(Priors,d_bins)
        #if np.all(Combined_bins!=Prior_bins):
        #    raise RuntimeError("Something went wrong with the binning")
    
    Combined_PDF = PDFs[0]
    for i in range(1,N):
        Combined_PDF *= PDFs[i]
        
    # Check that the Combined PDF is not 0
    check_combined(Combined_PDF,PDFs,Combined_bins,HIST="True")
    
    if Prior is not None:
        """Combined_Prior = Prior[0]
        for i in range(1,N):
            Combined_Prior *= Priors[i]
        #Check that the Combined Prior is not 0
        check_combined(Combined_Prior,KDE=True)
        Combined_Prior /= np.sum(Combined_Prior*get_bins_volume(d_bins)) # normalised 
        """
        check_combined(Prior,KDE=True)
        Norm_Prior = Prior/np.sum(Prior*get_bins_volume(d_bins))
        if Dim==3:
            plot = plot_probability3D(Norm_Prior,Combined_bins,labels=["AB","AC","AD"],udm="\"",alpha=0.3)
            plot.savefig(str(savedir)+"/Normalised_prior.png")
        for i in range(N-1):
            Combined_PDF/=Norm_Prior                
    return Combined_PDF, Combined_bins


def Multiply_PDF_KDE(samples,npoints,ratio=2,savedir="."):
    """
    Main idea for the moment:
        - get KDEs bandwiths
        - train KDEs at set of points from MCMC (maybe subsample of it)
        - set a set of points (~bin position) common for all KDEs
        - multiply KDEs at these points positions
        - we now have 1 KDE and a set of point where it has been evalueated
        - normalised point wise
    """
    #from sklearn.neighbors import KernelDensity
    #from sklearn.model_selection import GridSearchCV
    
    maxs,mins = get_minmax(samples)
    N         = len(samples)
    

    # select ony 1/n points of the sample -> select last n points from the center of the sample
    sparse_sample = []
    for smpl_fi in samples:
        points_smpl   = int(len(smpl_fi[0])/ratio)
        sparse_sample.append(np.array(smpl_fi).T[-points_smpl:].T)

        
    npoints   = complex(0,npoints)
    xx, yy,zz = np.mgrid[mins[0]:maxs[0]:npoints, mins[1]:maxs[1]:npoints, mins[2]:maxs[2]:npoints]
    positions = np.vstack([xx.ravel(), yy.ravel(),zz.ravel()])

    values    = [np.vstack(smpl) for smpl in sparse_sample]
    
    #kernels   = [st.gaussian_kde(np.transpose(vl)) for vl in values]
    #pool      = Pool(processes=100)  
    #kernels   = [pool.map(st.gaussian_kde,np.transpose(vl)) for vl in values]
    #KDEs      = [np.reshape(krnl(np.transpose(positions)).T, xx.shape) for krnl in kernels]

    with Pool(processes=100) as pool:
        kernels   = [pool.starmap(st.gaussian_kde,((np.transpose(vl),),))[0] for vl in values]
    
    #kernels   = [st.gaussian_kde(np.transpose(vl)) for vl in values]
    #KDEs      = [np.reshape(krnl(np.transpose(positions)).T, xx.shape) for krnll in kernels]
    KDEs      = [np.reshape(krnl(positions).T, xx.shape) for krnl in kernels]

    with open(savedir/str("PDF_KDEs.pkl"),"wb") as f:
        pickle.dump(KDEs,f)
    # we then have to sample the space

    
    Combined_PDF = KDEs[0]
    for i in range(1,N):
        Combined_PDF *= KDEs[i]
    check_combined(Combined_PDF,KDE=True)
    # normalise point-wise: Sum over all points =1
    Combined_PDF /= np.sum(Combined_PDF)

    # Prior computed here in the same position as the combined post -> only once

    if Prior is not None:
        # TEST: only once
        Prior,Positions_KDE_Prior = KDE_prior([Prior],positions,shape=xx.shape)

        Combined_Prior = Prior
        
        # Check that the Combined Prior is not 0
        check_combined(Combined_Prior,HIST="False")
        #Combined_Prior /= np.sum(Combined_Prior*get_bins_volume(d_bins)) # normalised 
        Combined_Prior /= np.sum(Combined_Prior) # normalised 
        plot = plot_probability3D_KDE(Combined_Prior,Positions_KDE_Prior,labels=["AB","AC","AD"],udm="\"",alpha=0.3)
        plot.savefig(str(savedir)+"/Combined_prior.pdf")
        for i in range(N-1):
            Combined_PDF/=Combined_Prior                
    return Combined_PDF, positions




def KDE_prior_corr(Prior, NtotPrior,d_bins,shape=None):
    # Priors shape must be (point_i, dimensions), eg: (1000,3)
    # Prior corrected for the Lower resolution
    # return Priors evaluated with KDE at center of the bins
    
    # find the best bandwidth for each filter
    prior_bndws = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth': prior_bndws})
    # initial shape of samples: filter,dimension,step

    bandwidth = grid.fit(Prior)
    #bandwidth_LR = grid.fit(Prior_LR)
    # input of grid.fit must have shape: ( n*steps, n*dimensions)

    # define the KDE
    kde = KernelDensity(**bandwidth.best_params_, kernel='gaussian') 
    #kde_LR = KernelDensity(**bandwidth_LR.best_params_, kernel='gaussian') 
    # Fit the KDE to the samples
    KDE = kde.fit(Prior)
    #KDE_LR = kde_LR.fit(Prior_LR)

    # we then have to sample the space
    if shape is None:
        positions = KDE_position_from_bin(d_bins) #position.shape = ((n_bins-1)**Dim,Dim) 
        shape     = [len(b)-1 for b in d_bins]
    else:
        positions = d_bins
    # NOTE 14th Nov '22:
    # the prior is computed by summing the high res prior to the prior obtained
    # from the low resolution 
    # They must have the same normalisation. First normalised "pointwise" both of them
    # then the high res have to be adapted to the approximately correct normalisation of the
    # low res prior. Then it have to be shifted st its average is 0
    # and finally it can be shifted to the real average given by the low res one
    Prior_val = KDE.score_samples(positions).reshape(*shape)
    #Prior_val_LR = KDE_LR.score_samples(positions).reshape(*shape)
    #shifted_PHR  = Prior_val - np.mean(Prior_val)
    #Prior        = shifted_PHR  + Prior_val_LR
    Prior = Prior_val* len(Prior)/NtotPrior
    return Prior, positions

def KDE_prior(Priors,d_bins,shape=None):
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
    if shape is None:
        positions = KDE_position_from_bin(d_bins) #position.shape = ((n_bins-1)**Dim,Dim) 
        shape     = [len(b)-1 for b in d_bins]
    else:
        positions = d_bins
        
    Priors = [kde_i.score_samples(positions).reshape(*shape) for kde_i in KDEs]
    
    # NOTE: by construction the KDE is not normalised
    # we could "normalise" it pointwise, but it is actually not necessary, 
    # as it is only shifted by a factor, and the Combined_PDF has to be normalised
    # anyway after the multiplication
    # but this is blocking the test so
    #Priors = [Prior_i/np.sum(Prior_i*get_bin_volume(d_bins)) for Prior_i in Priors] # "normalised" 
    return Priors, positions



"""
def wip_Multiply_PDF_HIST(samples,nbins,Prior=None,savedir="."):
    # Using histograms to obtain the binned density for each datasets
    N = len(samples)
    Dim = len(samples[0])
    #First find the same bins for each dataset in all parametric dimensions
    mins,maxs  = get_minmax(samples)#[DfABmax,DfACmax,DfBCmax],[DfABmin,DfACmin,DfBCmin]
    len_bins = (maxs-mins)/nbins
    d_bins = np.transpose([mins +i*len_bins for i in range(nbins+1)])
    
    # Then we make the histogram    
    PDFs,Combined_bins = hists_given_bins(samples,d_bins)
    
    if Prior is not None: # and Prior_LR is not None
        if Prior.uniform:
            # then it can be computed numerically (and it's a constant, so it doesn't really matter):
            Prior = np.ones_like(d_bins)/(nbins*Dim)
        else:
            Prior,_ = wip_KDE_prior_corr(Prior=Prior,d_bins=d_bins)
            
    Combined_PDF = PDFs[0]
    for i in range(1,N):
        Combined_PDF *= PDFs[i]
    
    # Check that the Combined PDF is not 0
    check_combined(Combined_PDF,PDFs,Combined_bins,HIST="True")
    
    if Prior is not None:
        check_combined(Prior,KDE=True)
        Norm_Prior = Prior/np.sum(Prior*get_bins_volume(d_bins))
        if Dim==3:
            plot = plot_probability3D(Norm_Prior,Combined_bins,labels=["AB","AC","AD"],udm="\"",alpha=0.3)
            plot.savefig(str(savedir)+"/Normalised_prior.png")
        for i in range(N-1):
            Combined_PDF/=Norm_Prior                
    return Combined_PDF, Combined_bins

def wip_KDE_prior_corr(Prior,d_bins,shape=None):
    # Prior is given as class and is considered Uniform
    # return Priors evaluated with KDE at center of the bins
    
    # find the best bandwidth for each filter
    prior_bndws = 10 ** np.linspace(-1, 1, 100)
    grid = GridSearchCV(KernelDensity(kernel='gaussian'),{'bandwidth': prior_bndws})
    # initial shape of samples: filter,dimension,step
    prior_sample = Prior.get_sample()
    bandwidth = grid.fit(prior_sample)
    # define the KDE
    kde = KernelDensity(**bandwidth.best_params_, kernel='gaussian') 
    
    # Fit the KDE to the samples
    KDE = kde.fit(prior_sample)
    
    # we then have to sample the space
    if shape is None:
        positions = KDE_position_from_bin(d_bins) #position.shape = ((n_bins-1)**Dim,Dim) 
        shape     = [len(b)-1 for b in d_bins]
    else:
        positions = d_bins

    Prior_val = KDE.score_samples(positions).reshape(*shape)
    Prior     = Prior_val* len(prior_sample)/Prior.NtotPrior
    return Prior, positions
"""

# General function to multiply PDF  (can be applied to priors as well - and will be)
"""def Multiply_PDF(samples,nbins=None,KDE=False,Priors=None,savedir="."):
    # samples.shape= (n_datasets,n_mcmc_points,n_dims)
    # same for Priors
    
    if KDE is False and nbins is None:
        raise RuntimeError("Nbins are required for the histogram sampling")

    if KDE:
        Combined_PDF,Positions = Multiply_PDF_KDE(samples,Priors,savedir)
        return Combined_PDF,Positions
    else:
        Combined_PDF, Combined_bins = Multiply_PDF_HIST(samples,nbins,Priors,savedir)
        return Combined_PDF, Combined_bins 
    
def get_minmax(samples,dim=None):
    # We get the min and maximum for each dimension relatively to all the datasets togheter
    if dim is not None:
        return np.min([min(par_i) for par_i in np.array(samples,dtype=object)[:,dim]]),\
            np.max([max(par_i) for par_i in np.array(samples,dtype=object)[:,dim]])
    else:
        dim = len(samples[0])
        min_i,max_i = [],[]
        for d in range(dim):
            min_i.append(np.min([min(par_i) for par_i in np.array(samples,dtype=object)[:,d]]))
            max_i.append(np.max([max(par_i) for par_i in np.array(samples,dtype=object)[:,d]]))
        return np.array(min_i),np.array(max_i)
"""

