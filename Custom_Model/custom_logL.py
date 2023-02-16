#!/usr/bin/env python
# coding: utf-8

import numpy as np

from Data.conversion import qphi_from_e1e2
from Data.input_data import init_kwrg_likelihood


class logL_ellipticity_aligned(object):

    def __init__(self, SUB, phi_ll=None):
        self.SUB    = SUB
        self.phi_ll = phi_ll
        pass
    
    def __call__(self, kwargs_lens=None, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)

    def logL_ellipt(self,phi,phi_ll,sigma_phi=4.5,bound=None):
        # give the log_L of a clipped normal likel.
        # all values are given in deg.
        logL_ell = -(phi-phi_ll)**2/sigma_phi**2/2 
        if bound:
            if abs(phi-phi_ll)>=bound:
                return  -10.**5 #return -np.inf
            else:
                return logL_ell
        else:
            return logL_ell

    
    def logL_addition(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        """
        a definition taking as arguments 
        (kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        if self.SUB:
            phi_ll = self.phi_ll
        else:
            e1_ll,e2_ll = kwargs_lens_light[0]["e1"],kwargs_lens_light[0]["e2"]
            q_ll,phi_ll = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
        
        e1,e2 = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
        
        logL  = self.logL_ellipt(phi,phi_ll)
        return logL

class logL_ellipticity_qphi(object):

    def __init__(self, SUB, phi_ll=None, q_ll=None):
        self.SUB    = SUB
        self.phi_ll = phi_ll
        self.q_ll   = q_ll
        pass
    
    def __call__(self, kwargs_lens=None, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)

    def logL_ellipt_phi(self,phi,phi_ll,sigma_phi=4.5,bound=None):
        # give the log_L of a clipped normal likel.
        # all values are given in deg.
        logL_ell = -(phi-phi_ll)**2/sigma_phi**2/2 
        if bound:
            if abs(phi-phi_ll)>=bound:
                return  -10.**5 #return -np.inf
            else:
                return logL_ell
        else:
            return logL_ell

    def logL_ellipt_q(self,q,q_ll,sigma_q=.1):
        # follow Note 27th June
        diff = q-(q_ll-0.1)
        if diff>=0:
            return 0
        else:
            return -diff**2/sigma_q**2/2
        

    def logL_addition(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        """
        a definition taking as arguments 
        (kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        if self.SUB:
            phi_ll = self.phi_ll
            q_ll   = self.q_ll
        else:
            e1_ll,e2_ll = kwargs_lens_light[0]["e1"],kwargs_lens_light[0]["e2"]
            q_ll,phi_ll = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
        
        e1,e2 = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
        
        logL  =  self.logL_ellipt_phi(phi,phi_ll)
        logL  += self.logL_ellipt_q(q,q_ll)
        return logL


#################### TEST #######################
# MOD_LLFR
from lenstronomy.LensModel.lens_model import LensModel
from Posterior_analysis.mag_remastered import get_mag_ratio

class logL_combined(object):
    def __init__(self, Rmag_ABC,sig_Rmag_ABC,lens_model,SUB, phi_ll=None, q_ll=None):
        # consider at same time ellipticity and magratio punishing terms:
        ######
        # Mag Ratio
        # note: only A,B and C are measured in the lcs, not D
        # note: input are RATIOS of magnification wrt to image A
        # note: the order is B/A,C/A
        if len(Rmag_ABC)>2:
            raise RuntimeError("More mag ratios then expected (2): ",len(Rmag_ABC))
        self.Rmag_ABC     = Rmag_ABC
        self.sig_Rmag_ABC = sig_Rmag_ABC
        self.lens_model   = lens_model
        ########
        # Ellipticity
        #
        self.SUB    = SUB
        self.phi_ll = phi_ll
        self.q_ll   = q_ll
        pass
    
    def __call__(self, kwargs_lens=None, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)

    def logL_ellipt_phi(self,phi,phi_ll,sigma_phi=4.5,bound=None):
        # give the log_L of a clipped normal likel.
        # all values are given in deg.
        logL_ell = -(phi-phi_ll)**2/sigma_phi**2/2 
        if bound:
            if abs(phi-phi_ll)>=bound:
                return  -10.**5 #return -np.inf
            else:
                return logL_ell
        else:
            return logL_ell

    def logL_ellipt_q(self,q,q_ll,sigma_q=.1):
        # follow Note 27th June
        diff = q-(q_ll-0.1)
        if diff>=0:
            return 0
        else:
            return -diff**2/sigma_q**2/2
    
    def logL_Rmag(Rmag_model,Rmag_obs,sig_Rmag_obs):
        Diff = np.array(Rmag_model)-np.array(Rmag_obs)
        sig2 = sig_Rmag_obs**2
        log  = -Diff**2/(2*sig2)
        return log
    
    def logL_addition(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        """
        a definition taking as arguments 
        (kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)
        and returns a logL (punishing) value.
        """
        logL = 0
        # Gaussian prior on the mag fraction obtained from the combined result 
        # of a specific lcs analysis            
        # abs: we ignore parity
        Rmag_model     = np.abs(get_mag_ratio(lens_model,kwargs_lens,kwargs_ps) )
        # only consider B/A, C/A:
        Rmag_model_ABC = Rmag_model[:-1] 
        logL  +=  self.logL_Rmag(Rmag_model_ABC,Rmag_obs=self.Rmag_ABC,sig_Rmag_obs=self.sig_Rmag_ABC)
        
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        if self.SUB:
            phi_ll = self.phi_ll
            q_ll   = self.q_ll
        else:
            e1_ll,e2_ll = kwargs_lens_light[0]["e1"],kwargs_lens_light[0]["e2"]
            q_ll,phi_ll = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
        
        e1,e2 = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
        
        logL  +=  self.logL_ellipt_phi(phi,phi_ll)
        logL  += self.logL_ellipt_q(q,q_ll)
        return logL


def init_kwrg_custom_likelihood(setting,mask=None,custom="qphi"):
    kwargs_likelihood = init_kwrg_likelihood(setting,mask)
    phi_ll = setting.phi_ll if setting.sub else None
    q_ll   = setting.q_ll   if setting.sub else None
    if custom=="aligned":
        kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_aligned(SUB=setting.sub,phi_ll=phi_ll)
    elif custom=="qphi":
        kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_qphi(SUB=setting.sub,phi_ll=phi_ll,q_ll=q_ll)
    elif custom=="combined":
        #pragma=no cover
        raise RuntimeError("the logL_ellipticiy_combined has to be implemented")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_combined()
    else:
        raise RuntimeError("The 'custom' string parameter must be :'aligned','qphi' or 'combined'(WIP). Not "+str(custom))
    return kwargs_likelihood