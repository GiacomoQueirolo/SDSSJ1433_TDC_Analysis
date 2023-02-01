#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from util.conversion_SDSS import qphi_from_e1e2
import numpy as np

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


# In[ ]:


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
        if diff<0:
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

