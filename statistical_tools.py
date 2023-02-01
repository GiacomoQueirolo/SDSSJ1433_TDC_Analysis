#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from tools import ProgressBar
# My sampling procedure

def inside_bin_i(pos,i,bins): 
    pos = posT(pos) #shape = n_dim, 1 
    inside=0
    for dim_i in range(len(pos)):
        if inside_1Dbin(x=pos[dim_i][0],bins_1D=bins[dim_i],index_bin=i):
            inside+=1
    return inside==len(pos)

def inside_1Dbin(x,index_bin,bins_1D):
    return x>=bins_1D[index_bin] and x<bins_1D[index_bin+1]
    
def posT(pos):
    # such that the shape of pos is now the same as the bin one
    if len(np.shape(pos))==2:
        return pos
    elif len(np.shape(pos))!=0:
        return [[p] for p in pos]
    else:
         return [[pos]]

def get_val_at_pos(pos,values,bins):
    coords = []
    """
    pos_T = posT(pos)
    for dim_i in range(len(pos_T)):        
        for bin_i in range(len(values[dim_i])):
            if inside_bin_i(pos_T,bin_i,bins):
                coords.append(bin_i)
                break
    """
    for dim_i in range(len(pos)):
        for index_bin in range(len(bins[dim_i])-1):
            if inside_1Dbin(pos[dim_i],index_bin,bins[dim_i]) :
                coords.append(index_bin)
                break    
    if len(coords)!=len(np.shape(values)):
        raise Exception("Position given not found in the binned region")
    return np.array(values).item(*coords)

def get_prob_at_pos(pos,values,bins):
    try:
        prob = get_val_at_pos(pos,values,bins)
    except Exception:
        prob = 0
    if prob>=0:
        return prob
    else:
        raise ValueError("Negative probability!")

def uniform_prior(prob):
    prior = np.ones_like(prob)/np.sum(np.ones_like(prob))
    return prior
def propose_mcmc_step(pos,mcmc_sigma,mcmc_range,n_prop_max=300):
    #most brute force way
    n_proposals = 0
    while True:
        n_proposals+=1
        new_pos = np.random.normal(pos,mcmc_sigma)
        if  pos_in_range(new_pos,mcmc_range):
            break
        elif n_proposals>n_prop_max:
            raise RuntimeError("Couldn't find a new mcmc position - ajust sigma or range")
    return new_pos

def pos_in_range(pos,mcmc_range):
    if np.shape(pos)==():
        pos = [pos]
        mcmc_range = [mcmc_range]
    for dim_i in range(len(pos)):
        if inside_1Dbin(pos[dim_i],0,mcmc_range[dim_i]):
            return False
    return True

def sampler(init_pos,prob,bins,mcmc_sigma,mcmc_steps,mcmc_range=None,prior="unif",text_mcmc="MCMC sampler"): # now generalised at n dim
    if prior =="unif":
        prior=uniform_prior(prob)
    if mcmc_range==None:
        mcmc_range=[]
        for dim_i in range(len(init_pos)): 
            mcmc_range.append( [np.min(bins[dim_i]),np.max(bins[dim_i])])
    pos       = init_pos
    like_pos  = get_prob_at_pos(init_pos,prob ,bins)
    prior_pos = get_prob_at_pos(init_pos,prior,bins)
    prob_pos  = like_pos*prior_pos
    mcmc_chain = [pos]
    mcmc_likelihood = [like_pos]
    n_step = 0
    while n_step<mcmc_steps:
        pos_proposed   = propose_mcmc_step(pos,mcmc_sigma,mcmc_range) #np.random.normal(pos,mcmc_sigma) # set a range
        like_proposed  = get_prob_at_pos(pos_proposed,prob ,bins)
        prior_proposed = get_prob_at_pos(pos_proposed,prior,bins)
        
        prob_proposed  = like_proposed*prior_proposed
        if prob_pos==0.0:   
            prob_pos=1e-200        
        p_accept = prob_proposed/prob_pos
            
        accept  =  p_accept > np.random.random()
        if accept:
            pos = pos_proposed
            like_pos = like_proposed
            mcmc_chain.append(np.array(pos).tolist())
            mcmc_likelihood.append(like_pos)
        n_step+=1
        ProgressBar(n_step,mcmc_steps,text_mcmc)
    chain = [np.array(mcmc_chain),np.array(mcmc_likelihood)]
    return chain


# In[ ]:


def get_bins_volume(bins):    
    vol = []
    for dim_i in range(len(bins)):
        vol_i = []
        for bin_i in range(len(bins[dim_i])-1):
            vol_i.append(bins[dim_i][bin_i+1]-bins[dim_i][bin_i])
        vol.append(vol_i)
    if len(vol)!=3:
        raise RuntimeError("Implemented only for 3D")     ### To resolve
    vol_grid=[]
    for bin_i in range(len(vol[0])):
        vol_bin_i = []
        for bin_j in range(len(vol[1])):
            vol_bin_j = []
            for bin_k in range(len(vol[2])):
                vol_bin_k = vol[0][bin_i]*vol[1][bin_j]*vol[2][bin_k]
                vol_bin_j.append(vol_bin_k)
            vol_bin_i.append(vol_bin_j)
        vol_grid.append(vol_bin_i)
    return np.array(vol_grid)
   
#def vector_mult(Vects,Index):
#    return np.prod([Vects[i][index] for i,index in enumerate(Index)]) 

def marginalise_prob(prob,bins):
    #marginalise over the other dimensions
    marg_prob = []
    vol_bins = get_bins_volume(bins)
    all_axis = np.arange(len(np.shape(prob)))
    for i in range(len(all_axis)):
        ax = tuple(np.delete(all_axis,i).tolist())
        marg_prob_i = np.sum(prob*vol_bins,axis=ax)
        marg_prob.append(marg_prob_i)
    
    for i in range(len(marg_prob)):
        if abs(np.sum(marg_prob[i])-1)>1e-5:
            # This method won't work if the prob is not normalised (quantiles!) 
            raise ValueError("Not normalised marginalised prob: sum=",np.sum(marg_prob[i]))
    return marg_prob
def mean_bins(bins,bin_i):
    # return mean of two consecutive bins given the first @ bin_i
    return 0.5*(bins[bin_i]+bins[bin_i+1])

def estimate_median(prob_marg,bins):
    median = []
    for dim_i in range(len(prob_marg)):
        median_i = None
        qnt = 0.
        for bin_i in range(len(prob_marg[dim_i])):
            qnt += prob_marg[dim_i][bin_i]
            if qnt>=0.5:
                median_i = mean_bins(bins[dim_i],bin_i)
                break
        median.append(median_i)
    return median

def estimate_sigma(prob_marg,bins,median=None,averaged=True):
    if median is None:
        median = estimate_median(prob_marg,bins)
    sigma = []
    for dim_i in range(len(prob_marg)):
        sigma_low  = None
        sigma_high = None
        qnt = 0.
        for bin_i in range(len(prob_marg[dim_i])):
            qnt += prob_marg[dim_i][bin_i]
            if qnt>=0.16 and sigma_low is None:
                sigma_low = abs(median[dim_i] -mean_bins(bins[dim_i],bin_i))
            elif qnt>=0.84 and sigma_high is None:
                sigma_high = abs(mean_bins(bins[dim_i],bin_i)- median[dim_i])
        if averaged:
            sigma_i = .5*(sigma_high+sigma_low)
            if sigma_i<1e-20:
                sigma_i = 1e-20
        else:
            sigma_i = [sigma_low,sigma_high]
        sigma.append(sigma_i)
    return sigma

def estimate_for_mcmc(prob,bins):
    prob_marg  = marginalise_prob(prob,bins)
    median     = estimate_median(prob_marg,bins)
    sigma      = estimate_sigma(prob_marg,bins,median=median)
    return median,sigma

# for PH0 calculation
def quantiles(prob,sampling_prob,q=[0.16,.5,0.84],return_quantiles=False):
    qnt    = []
    integr = 0
    for i in range(len(prob)):
        integr+=prob[i]
        for qi in q:
            if integr-prob[i]<qi and integr>=qi:
                qnt.append(sampling_prob[i])
    if return_quantiles:
        return qnt
    else:
        res     = qnt[1]
        err_min = qnt[1]-qnt[0]
        err_max = qnt[2]-qnt[1]
        return (res,err_min,err_max)

