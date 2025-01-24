# To use in conjunction with Prior.prior_likelihood_KDE

# We want it to have a flat prior for the ellipticity parameters q,phi
# which means a not flat prior for e1,e2 BUT flat for all remaining params
# so we sample the KDE for e1,e2, but make it so that the other params are only valued
# within the ranges 

import numpy as np
from sklearn.neighbors import KernelDensity

from Utils.tools import *
from Data.conversion import e1e2_from_qphi

class KDE_prior():
    def __init__(self,KDE):
        self.kde=KDE
    def likelihood(self,pos):
        logL = self.kde.score_samples(pos)# this is in log density
        return float(np.exp(logL))
    
class UNIF_prior():
    def __init__(self,bound_min,bound_max):
        self.bound_min = bound_min
        self.bound_max = bound_max
    def likelihood(self,pos):
        likl = 1 
        if pos < self.bound_min or pos >self.bound_max:
            likl = 0
        return likl

"""# Sample from a uniform distribution
uniform_samples = np.random.uniform(0, 1, 60000)

# Transform the samples using the inverse CDF
samples = inverse_cdf(uniform_samples,x0=0,sig=sig_min)

# Plot the original function and the sampled distribution
plt.plot(x_lnsp, semigauss_norm_arr(x_lnsp,x0=0,sig=sig_min), label='Original Function')
plt.hist(samples, bins=50, density=True, alpha=0.5, label='Sampled Distribution')
plt.legend()
plt.show()"""
class inefficient_semigaussian_q_sampling():
    def __init__(self,savepath=None,n_x=60000):
        self.n_x      = n_x
        self.q_lnsp   = np.linspace(0, 1, n_x)
        self.savepath = savepath
    def semigauss(self,x,x0,sig):
        if x<x0:
            return np.exp(-(x-x0)**2/(2*sig**2))
        else:
            return 1
    def semigauss_arr(self,x,x0,sig):
        return np.array([self.semigauss(xi,x0,sig) for xi in x])

    def semigauss_norm_arr(self,x,x0,sig):
        y = self.semigauss_arr(x,x0,sig) 
        norm = np.trapz(self.semigauss_arr(self.q_lnsp,x0,sig),self.q_lnsp)
        return y/norm
    # Calculate the cumulative distribution function (CDF)
    def cdf(self,q,q0,sig):
        #return np.trapz([semigauss_norm_arr(xi,x0=x0,sig=sig) for xi in np.linspace(-np.inf, x, 1000)], dx=0.001)
        xx  = np.linspace(0, q, self.n_x)
        smg = self.semigauss_norm_arr(xx,x0=q0,sig=sig)
        return np.trapz(smg,xx)

    # Define the inverse of the CDF
    def inverse_cdf(self,y,q0,sig):
        cdf_comp = [self.cdf(qi,q0,sig) for qi in self.q_lnsp]
        return  np.array([float(self.q_lnsp[np.where(np.abs(cdf_comp-i)==min(np.abs(cdf_comp-i)))]) for i in y])
    
    def get_samples(self,q0,sig,unif_samples=None,savepath=None):
        if savepath is None:
            savepath = self.savepath
            try:
                with open(savepath,"rb") as f:
                    _samples = pickle.load(f)
                if len(_samples)!=self.n_x:
                    raise RuntimeError("Lenght of previously computed prior is different, recomputing")
                else:
                    self.samples = _samples
                return self.samples
            except:
                pass
        if hasattr(self,"samples"):
            return self.samples
        if not unif_samples:
            unif_samples = np.random.uniform(0, 1,self.n_x)
        self.samples = self.inverse_cdf(unif_samples,q0=q0,sig=sig)
        if savepath is not None:
            with open(savepath,"wb") as f:
                pickle.dump(self.samples,f) 
        return self.samples


class Semigauss_qphi_KDE_prior():
    # this is SPECIFIC only for the MAIN lens with q,phi as "semigaussian"
    # all the remaining are considered flat
    # see notes 14th Nov. '23
    def __init__(self,kwargs_min,kwargs_max,kwargs_fix,kwargs_qphi,nsample=int(5e4)):
        #self.kwargs_lens = kwargs_main_lens
        # shape : kw_init,kw_sig,kw_fx,kw_min,kw_max
        self.kwargs_min,self.kwargs_max = kwargs_min,kwargs_max
        self.kwargs_fix                 = kwargs_fix
        self.kwargs_qphi = kwargs_qphi        
        self.nsample     = nsample
        self.savepath_qsamples =  None
        if "savepath_qsamples" in kwargs_qphi:
            self.savepath_qsamples = kwargs_qphi["savepath_qsamples"]

        self.kde_e1e2_prior = self.get_e1e2_KDE(kwargs_qphi)

    
    def _get_qphi_semigaussian_prior(self,kwargs_qphi):
        q0,sig_q0   = kwargs_qphi["q"]        
        q_prior_cl  = inefficient_semigaussian_q_sampling(savepath=self.savepath_qsamples,n_x=self.nsample)
        q_prior     = q_prior_cl.get_samples(q0=q0,sig=sig_q0)
        phi0,sig_phi0 = kwargs_qphi["phi"]
        phi_prior     = np.random.normal(phi0,sig_phi0,self.nsample)
        return q_prior,phi_prior

    def get_e1e2_KDE(self,kwargs_qphi):
        if hasattr(self,"kde_e1e2_prior"):
            return self.kde_e1e2_prior
        q_prior,phi_prior   = self._get_qphi_semigaussian_prior(kwargs_qphi)
        e1_pr,e2_pr         = e1e2_from_qphi(q_prior,phi_prior,deg=self.kwargs_qphi["phi_degree"])
        self.kde_e1e2_prior = get_kde_prior(np.transpose([e1_pr,e2_pr]))
        return self.kde_e1e2_prior
    """    
    def _index_e1e2(self):
        keys = list(self.kwargs_max.keys())
        return keys.index("e1"),keys.index("e2")
    def likelihood(self,pos):
        # this pos is a vector 
        if len(pos)!=len(self.kwargs_max):
            raise RuntimeError("pos has to have the same lenght as the main lens parameters")
        likl = 1 
        for i in range(len(pos)):
            if self.kwargs_max.keys()[i]!="e1" and self.kwargs_max.keys()[i]!="e2":
                if pos < self.kwargs_min[i] or pos >self.kwargs_max[i]:
                    return 0 
        index_e1e2 = self._index_e1e2()
        e1e2_pos = [pos[index_e1e2[0]],pos[index_e1e2[1]]]
        lkl_qphi = self.kde_e1e2_prior.likelihood(e1e2_pos)
        likl*=lkl_qphi
        return likl
    """
    def likelihood(self,kw_pos):
        # this pos is a kwargs containing ALL params value of the lens parameters (all lenses)

        likl = 1 
        
        for i in range(len(kw_pos)):
            for prm in kw_pos[i].keys():
                if  "e1"  in self.kwargs_max[i].keys() and "e2"  in self.kwargs_max[i].keys():
                    index_ml = i #index of main lens model
                    if  prm=="e1" or prm=="e2":
                        continue
                if prm in self.kwargs_fix[i].keys():
                    continue
                elif kw_pos[i][prm] < self.kwargs_min[i][prm] or kw_pos[i][prm]>self.kwargs_max[i][prm]:
                    #if it's out of bounds
                    return 0               
        e1e2_pos = np.array([kw_pos[index_ml]["e1"],kw_pos[index_ml]["e2"]]).reshape(1,-1)
        lkl_qphi = self.kde_e1e2_prior.likelihood(e1e2_pos)
        likl    *= lkl_qphi
        return likl
    





def get_kde_prior(Prior):
    #bandwidth = grid.fit(Prior)
    #kde = KernelDensity(**bandwidth.best_params_, kernel='gaussian') 
    # hiperparams tuned "by hand" in remade_ellipt_recentered
    kde = KernelDensity(bandwidth=0.03, kernel='tophat')  
    kde.fit(Prior)
    return KDE_prior(kde)



class Precomputed_KDE():
    def __init__(self,setting,kwrg_qphi={"type":"Semigaussian"},nsample=int(5e4)):
        self.setting         = setting
        self.kwargs_params   = get_kwargs_params(setting)
        self._nsample        = nsample
        self.kwrg_qphi       = kwrg_qphi
        if self.kwrg_qphi is None:
            self.lens_list       = self._get_unif_prior_list(self.kwargs_params["lens_model"])
        elif self.kwrg_qphi["type"]=="Gaussian":
            raise NotImplementedError("Not yet implemented the pure gaussian")
        elif self.kwrg_qphi["type"]=="Semigaussian":
            if "q" not in kwrg_qphi.keys():
                self.kwrg_qphi["q"]  = setting.pll["q"]
            if "phi" not in kwrg_qphi.keys():
                self.kwrg_qphi["phi"]  = setting.pll["phi"]
                self.kwrg_qphi["phi_degree"] = True
            self.lens_list      = self._get_semigaussian_qphi_prior_list(self.kwargs_params["lens_model"])
        self.lens_light_list = self._get_unif_prior_list(self.kwargs_params["lens_light_model"])
        self.ps_list         = self._get_unif_prior_list(self.kwargs_params["point_source_model"])
        self.source_list     = self._get_unif_prior_list(getattr(self.kwargs_params,"source_model",None))
    def _get_unif_prior_list(self,kwrg_prm):
        if kwrg_prm is None:
            return 0
        prior_list = []
        kw_in,kw_sg,kw_fx,kw_min,kw_max= kwrg_prm
        for i in range(len(kw_min)):
            for prm_nm in kw_max[i].keys():
                prior_list.append(UNIF_prior(bound_max=kw_max[i][prm_nm],bound_min=kw_min[i][prm_nm]))
        return prior_list
    def _get_semigaussian_qphi_prior_list(self,kwrg_prm):
        if kwrg_prm is None:
            return 0 
        kw_in,kw_sg,kw_fx,kw_min,kw_max= kwrg_prm
        return Semigauss_qphi_KDE_prior(kwargs_min=kw_min,kwargs_max=kw_max,kwargs_fix=kw_fx,\
                                 kwargs_qphi=self.kwrg_qphi,nsample=self._nsample)
        """
        # this has to be completely different, 
        # ie has to take a vector for the lens params and return the likelihood (not a list)
        for i_mod in range(len(kw_min)):
            if i_mod==0: #main lens
                prior_list.append(Semigauss_qphi_KDE_prior(kwargs_max=kw_max[i_mod],kwargs_min=kw_min[i_mod],\
                                                           kwargs_qphi=self.kwrg_qphi,nsample=self._nsample))
            else:
                for prm_nm in kw_max[i_mod].keys():           
                    prior_list.append(UNIF_prior(bound_max=kw_max[i_mod][prm_nm],bound_min=kw_min[i_mod][prm_nm]))
        return prior_list
        """
