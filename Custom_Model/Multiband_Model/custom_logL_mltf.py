from Data.conversion import qphi_from_e1e2
from Data.Multiband_Data.input_data_mltf import init_kwrg_likelihood_mltf
from Custom_Model.custom_logL import logL_ellipt_q,logL_ellipt_phi,logL_center_lens


class logL_Prior_Lens_Light_Mltf(object):

    def __init__(self, multifilter_setting):
        self.pll    = multifilter_setting.pll
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, \
                 kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens)

    def logL_addition(self, kwargs_lens):
        """
        a definition taking as arguments kwargs_lens
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        # moreover we consider a free center of the lens light and mass, 
        # but normally correlated
        
        e1,e2   = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  =  logL_ellipt_phi(phi=phi,phi_ll=self.pll["phi"][0])
        logL_add  += logL_ellipt_q(q=q,q_ll=self.pll["q"][0],sigma_q=self.pll["q"][1]*3)
        logL_add  += logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        return logL_add


#MOD_LLCNT: 
class logL_PLL_cnt_free_Mltf(object):
    # similar to PLL, but if the lens light is modelled (for infrareds),
    # it tight it's center and ellipt to PLL as well
    def __init__(self, multifilter_setting):
        self.pll  = multifilter_setting.pll
        self.subs = [sett.sub for sett in multifilter_setting.settings] 
        self.qll_sigma = getattr(multifilter_setting,"qll_sigma",0.17)
        
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
        
        e1,e2   = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  =  logL_ellipt_phi(phi=phi,phi_ll=self.pll["phi"][0])    
        logL_add  += logL_ellipt_q(q=q,q_ll=self.pll["q"][0],sigma_q=self.qll_sigma)
        logL_add  += logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        # check if main lens is done
        llm = 0
        for i in range(len(self.subs)):
            if self.subs[i] is False:
                e1_mll,e2_mll = kwargs_lens_light[llm]["e1"],kwargs_lens_light[llm]["e2"]
                q_mll,phi_mll = qphi_from_e1e2(e1_mll,e2_mll,ret_deg=True)
                x_mll,y_mll   = kwargs_lens_light[llm]["center_x"],kwargs_lens_light[llm]["center_y"]
                logL_add     +=  logL_ellipt_phi(phi=phi_mll,phi_ll=self.pll["phi"][0])    
                logL_add     += logL_ellipt_q(q=q_mll,q_ll=self.pll["q"][0],sigma_q=self.qll_sigma)
                logL_add     += logL_center_lens(xll=xll,yll=yll,xm=x_mll,ym=y_mll)
                llm+=3
            else:
                llm+=2
        return logL_add

class logL_Prior_Center_LL_Mltf(object):

    def __init__(self, multifilter_setting):
        self.pll    = multifilter_setting.pll
        
    def __call__(self, kwargs_lens, kwargs_source=None, kwargs_lens_light=None, \
                 kwargs_ps=None, kwargs_special=None, kwargs_extinction=None):
        return self.logL_addition(kwargs_lens)

    def logL_addition(self, kwargs_lens):
        """
        a definition taking as arguments kwargs_lens
        and returns a logL (punishing) value.
        """
        # Gaussian prior on the ellipticity of the lens mass profile given the corresponding 
        # light profile
        #  if LL is subtracted, the phi_ll is given as a fixed value
        #  else it is computed from the corresponding lens light model
        # moreover we consider a free center of the lens light and mass, 
        # but normally correlated
        
        e1,e2   = kwargs_lens[0]["e1"],kwargs_lens[0]["e2"]
        q,phi   = qphi_from_e1e2(e1,e2,ret_deg=True)
        xm,ym   = kwargs_lens[0]["center_x"],kwargs_lens[0]["center_y"]
        xll,yll = self.pll["x"][0],self.pll["y"][0]
        
        logL_add  += logL_center_lens(xll=xll,yll=yll,xm=xm,ym=ym)
        return logL_add




def init_kwrg_custom_likelihood_mltf(multifilter_setting,mask=None,custom="qphi"):
    kwargs_likelihood = init_kwrg_likelihood_mltf(multifilter_setting,mask)
    
    if custom=="aligned":
        raise RuntimeError(f"pragma no cover:{custom}")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_aligned(setting=setting)
    elif custom=="centered":
        kwargs_likelihood["custom_logL_addition"] = logL_Prior_Center_LL_Mltf(setting=setting)

    elif custom=="qphi":
        raise RuntimeError(f"pragma no cover:{custom}")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_qphi(setting=setting)
    elif custom=="PLL":
        kwargs_likelihood["custom_logL_addition"] = logL_Prior_Lens_Light_Mltf(multifilter_setting=multifilter_setting)
    elif custom=="PLL_LLCNT":
        kwargs_likelihood["custom_logL_addition"] = logL_PLL_cnt_free_Mltf(multifilter_setting=multifilter_setting)
    elif custom=="combined":
        #pragma=no cover
        raise RuntimeError("the logL_ellipticiy_combined has to be implemented")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_combined()
    elif custom==None:
        return kwargs_likelihood
    else:
        raise RuntimeError("The 'custom' string parameter must be :'aligned','qphi' or 'combined'(WIP). Not "+str(custom))
    return kwargs_likelihood
