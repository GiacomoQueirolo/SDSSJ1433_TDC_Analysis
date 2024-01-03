#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from corner import quantile
from Utils.tools import ProgressBar
"""
def get_median_with_error(prob,bins,ret_str=True):
    median = estimate_median([prob],[bins])[0]
    sigmas = estimate_sigma([prob],[bins],median=[median],averaged=False)[0]
    rounding = int(1-np.log10(max(sigmas)))
    if ret_str:
        sig_up_rnd,sig_down_rnd = str(np.round(sigmas[0],rounding)),str(np.round(sigmas[1],rounding))
        if sig_up_rnd==sig_down_rnd:
            return str(np.round(median,rounding))+"$\pm$"+sig_up_rnd
        return str(np.round(median,rounding))+"$_{"+sig_down_rnd+"}^{"+sig_up_rnd+"}$"
    else:
        return median,sigmas
        """
def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return float('.'.join([i, (d+'0'*n)[:n]]))

def simp(val_raw,sig_raw1,sig_raw2,ignore_debug=False): 
    _, fact1 =simp_sig(sig_raw1)
    _, fact2 =simp_sig(sig_raw2)
    
    if val_raw ==0:
        val= "0.0"
    elif fact1 is None and fact2 is None :
        if abs(val_raw)>1000 or abs(val_raw)<0.001:
            val= "{:.2e}".format(val_raw) 
        else:
            fact=1000
            val= str(np.round(val_raw,fact))
    else:
        fact = max(fact1,fact2)
        val  = np.round(val_raw,fact)
        sig1 = np.round(sig_raw1,fact)
        sig2 = np.round(sig_raw2,fact)
        if abs(val_raw)>1000 or abs(val_raw)<0.001:
            power_v  = -int(('%E' % val_raw) [-3:])
            power_s1 = -int(('%E' % sig_raw1)[-3:])
            power_s2 = -int(('%E' % sig_raw2)[-3:])
            power_s  = max([power_s1,power_s2])
            dp       = power_s - power_v
            val      = str(np.round(val_raw*(10**power_v),dp))
            sig1     = str(np.round(sig_raw1*(10**(power_v)),dp))
            sig2     = str(np.round(sig_raw2*(10**(power_v)),dp))
            if sig1!=sig2:        
                return "$ "+val+"^{+"+sig1+"}_{-"+sig2+r"}\cdot 10^{"+str(-power_v)+r"}$"
            else:
                return "$ "+val+ r"\pm"+sig1+r"\cdot 10^{"+str(-power_v)+r"}$"
        else:
            val =str(val)
            sig1=str(sig1) 
            sig2=str(sig2)
    if val=="0.0" and val_raw!=0 or len(sig1)>len(val) or len(sig2)>len(val) :
        if not ignore_debug:
            raise RuntimeError(f"#DEBUG\nval {val}\nval_raw {val_raw}\nsig_raw1 {sig_raw1}\nsig_raw2 {sig_raw2}\nsig1 {sig1}\nsig2 {sig2}\n#DEBUG")
    if sig1!=sig2:
        return "$ "+val + "^{+"+sig1+"}"+"_{-"+sig2+"} $"
    else:
        return "$ "+val + "\\pm"+sig2+" $"
            

def simp_wo_sig(val_raw,rnd_vl=1):
    if abs(val_raw)>1000 or abs(val_raw)<0.001:
        power_v  = -int(('%E' % val_raw) [-3:])
        val      = str(np.round(val_raw*(10**power_v),1))
        return "$ "+val+ r"\cdot 10^{"+str(-power_v)+r"}$"
    else:
        return "$"+str(np.round(val_raw,rnd_vl))+"$"


def simp_sig(val):
    if val ==0:
        return ("0.0",None)
    else: 
        fact=1
        while truncate(val,fact)==0.0:
            fact+=1
        return (str(np.round(val,fact)),fact)


def easy_res(mcmc,quantiles=None):
    median = np.median(mcmc)
    if not quantiles:
        return {"median":median,"std":np.std(mcmc)}
    else:
        qnt = quantile(mcmc,q=quantiles)
        min_err = median - np.min(qnt)
        max_err = -median + np.max(qnt)
        return {"median":median,"quantiles":quantiles,"min_err":min_err,"max_err":max_err}
        
def easy_res_str(mcmc,quantiles:bool=True):
    if quantiles:
        quantiles = [0.16,0.84]
    res = easy_res(mcmc,quantiles=quantiles)
    if not quantiles:
        rounding = int(1-np.log10(res["std"]))
        return f"${np.round(res['median'],rounding)}\pm{np.round(res['std'],rounding)}$"
    else:
        #rounding = int(1-np.log10(np.max([res["min_err"],res["max_err"]])))
        #return f"${np.round(res['median'],rounding)}"+"_{-"+str(np.round(res['min_err'],rounding))+"}^{+"+str(np.round(res['max_err'],rounding))+"}$"
        return simp(res['median'],res['min_err'],res["max_err"])

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

