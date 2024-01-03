#!/usr/bin/env python
# coding: utf-8

import numpy as np

from Data.conversion import qphi_from_e1e2
from Data.input_data import init_kwrg_likelihood



class logL_ellipticity_aligned(object):

    def __init__(self,setting):
        self.setting = setting
        self.SUB     = setting.SUB
        self.phi_ll  = setting.phi_ll
    
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
        

def logL_ellipt_phi(phi,phi_ll,sigma_phi=4.5,bound=None):
        # give the log_L of a clipped normal likel.
        # all values are given in deg.
        logL_phi = -(phi-phi_ll)**2/sigma_phi**2/2 
        if bound:
            if abs(phi-phi_ll)>=bound:
                return  -10.**5 #return -np.inf
            else:
                return logL_phi
        else:
            return logL_phi

def logL_ellipt_q(q,q_ll,sigma_q=.1):
        # follow Note 27th June
        diff = q-(q_ll-0.1)
        if diff>=0:
            return 0
        else:
            return -diff**2/sigma_q**2/2      

#DEPRECATED
class logL_ellipticity_qphi(object):

    def __init__(self,setting):
        self.setting = setting
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens,kwargs_lens_light)

    def logL_addition(self, kwargs_lens, kwargs_lens_light=None):
        """
        a definition taking as arguments 
        (kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)
        and returns a logL (punishing) value.
        """
        raise RuntimeWarning("Deprecated: this introduced biased mass model likelihood depending on the single filter lens light distribution")
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        if self.setting.sub:
            try:
                phi_ll = self.setting.phi_ll
                q_ll   = self.setting.q_ll
            except AttributeError:
                phi_ll = self.setting.pll["phi"][0]
                q_ll   = self.setting.pll["q"][0]
        else:
            e1_ll,e2_ll = kwargs_lens_light[0]["e1"],kwargs_lens_light[0]["e2"]
            q_ll,phi_ll = qphi_from_e1e2(e1_ll,e2_ll,ret_deg=True)
        
        e1,e2 = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
        logL_add  =   0
        logL_add  += logL_ellipt_phi(phi,phi_ll)
        logL_add  += logL_ellipt_q(q,q_ll)
        return logL_add
    

#################### TEST #######################
# MOD_LLFR
"""
from lenstronomy.LensModel.lens_model import LensModel
from Data.mag_remastered import get_mag_ratio

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
    
    def __call__(self, kwargs_lens=None, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)

    
    def logL_Rmag(self,Rmag_model,Rmag_obs,sig_Rmag_obs):
        Diff = np.array(Rmag_model)-np.array(Rmag_obs)
        sig2 = np.array(sig_Rmag_obs)**2
        log  = -Diff**2/(2*sig2)
        return log
    
    def logL_addition(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
"""
#a definition taking as arguments 
#(kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps, kwargs_special, kwargs_extinction)
#and returns a logL (punishing) value.
"""
        logL = 0
        # Gaussian prior on the mag fraction obtained from the combined result 
        # of a specific lcs analysis            
        # abs: we ignore parity
        # to verify: image order
        im_order = image_order(kwargs_ps[0]["ra_image"],kwargs_ps[0]["dec_image"],verbose=False)
        if im_order[1:3].tolist()!=[1,2]:
            kwargs_ps[0]["ra_image"] = np.array(kwargs_ps[0]["ra_image"])[im_order]
            kwargs_ps[0]["dec_image"] = np.array(kwargs_ps[0]["dec_image"])[im_order]
            
        Rmag_model = np.abs(get_mag_ratio(self.lens_model,kwargs_lens,kwargs_ps) )

        # only consider B/A, C/A:
        Rmag_model_ABC = Rmag_model[:-1] 
        logL_ABC = self.logL_Rmag(Rmag_model=Rmag_model_ABC,Rmag_obs=self.Rmag_ABC,sig_Rmag_obs=self.sig_Rmag_ABC)
        logL += np.sum(logL_ABC)
        
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
        
        logL  += logL_ellipt_phi(phi,phi_ll)
        logL  += logL_ellipt_q(q,q_ll)
        return logL

"""
    
##### TEST2 -> reworked for PLL
# considering the center of the lens light and lens mass as not fixed togheter but correlated

def logL_center_lens(xll,yll,xm,ym,sigma_d=.4):
        # note: distance in arcsec
        dx2  = (xll-xm)**2
        dy2  = (yll-ym)**2
        diff = np.sqrt(dx2+dy2)
        return -diff**2/sigma_d**2/2
        
class logL_Prior_Lens_Light(object):

    def __init__(self, setting):
        self.pll    = setting.pll
        # MOD_STRCTQ
        self.qll_sigma = getattr(setting,"qll_sigma",0.17)
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, \
                 kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens)

    def logL_addition(self, kwargs_lens):
        """
        a definition taking as arguments kwargs_lens
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile prior
        
        e1,e2   = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  =   0
        logL_add  += logL_ellipt_phi(phi=phi,phi_ll=self.pll["phi"][0])
        logL_add  += logL_ellipt_q(q=q,q_ll=self.pll["q"][0],sigma_q=self.qll_sigma)#self.pll["q"][1]*3)
        logL_add  += logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        return logL_add


#MOD_LLCNT: DEPRECATED
class logL_PLL_cnt_free(object):
    # similar to PLL, but if the lens light is modelled (for infrareds),
    # it tight it's center and ellipt to PLL as well
    def __init__(self, setting):
        self.pll = setting.pll
        self.sub = setting.sub 
        self.qll_sigma = getattr(setting,"qll_sigma",0.17)
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, \
                 kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens,kwargs_lens_light=kwargs_lens_light)

    def logL_addition(self, kwargs_lens,kwargs_lens_light=None):
        """
        a definition taking as arguments kwargs_lens
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile 
        # and light profile of the main lens given the corresponding 
        # initial fixed values
        
        raise RuntimeWarning("Deprecated: no need to correlated the Sersic lens light with the isophote one, they are independent models.")
        e1,e2   = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  =  0
        logL_add  += logL_ellipt_phi(phi=phi,phi_ll=self.pll["phi"][0])    
        logL_add  += logL_ellipt_q(q=q,q_ll=self.pll["q"][0],sigma_q=self.qll_sigma)
        logL_add  += logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        # check if main lens is done
        if not self.sub:
            e1_mll,e2_mll = kwargs_lens_light[0]["e1"],kwargs_lens_light[0]["e2"]
            q_mll,phi_mll = qphi_from_e1e2(e1_mll,e2_mll,ret_deg=True)
            x_mll,y_mll   = kwargs_lens_light[0]["center_x"],kwargs_lens_light[0]["center_y"]
            logL_add     +=  logL_ellipt_phi(phi=phi_mll,phi_ll=self.pll["phi"][0])    
            logL_add     += logL_ellipt_q(q=q_mll,q_ll=self.pll["q"][0],sigma_q=self.qll_sigma)
            logL_add     += logL_center_lens(xll=xll,yll=yll,xm=x_mll,ym=y_mll)
        
        return logL_add

# only the center
class logL_Prior_Center_LL(object):

    def __init__(self, setting):
        self.pll = setting.pll
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, \
                 kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens)

    def logL_addition(self, kwargs_lens):
        """
        a definition taking as arguments kwargs_lens
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the center of the lens mass profile given the corresponding 
        # light profile prior

        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  =    0
        logL_add  +=  logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        return logL_add


def init_kwrg_custom_likelihood(setting,mask=None,custom="qphi"):
    kwargs_likelihood = init_kwrg_likelihood(setting,mask)
    
    if custom=="aligned":
        kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_aligned(setting=setting)
    elif custom=="centered":
        kwargs_likelihood["custom_logL_addition"] = logL_Prior_Center_LL(setting=setting)
    elif custom=="qphi":
        kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_qphi(setting=setting)
    elif custom=="PLL":
        kwargs_likelihood["custom_logL_addition"] = logL_Prior_Lens_Light(setting=setting)
    elif custom=="PLL_LLCNT":
        kwargs_likelihood["custom_logL_addition"] = logL_PLL_cnt_free(setting=setting)
    elif custom=="combined":
        #pragma=no cover
        raise RuntimeError("the logL_ellipticiy_combined has to be implemented")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_combined()
    elif custom==None:
        return kwargs_likelihood
    else:
        raise RuntimeError("The 'custom' string parameter must be :'aligned','qphi' or 'combined'(WIP). Not "+str(custom))
    return kwargs_likelihood

