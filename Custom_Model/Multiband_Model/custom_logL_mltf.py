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


def init_kwrg_custom_likelihood_mltf(multifilter_setting,mask=None,custom="qphi"):
    kwargs_likelihood = init_kwrg_likelihood_mltf(multifilter_setting,mask)
    
    if custom=="aligned":
        raise RuntimeError(f"pragma no cover:{custom}")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_aligned(setting=setting)
    elif custom=="qphi":
        raise RuntimeError(f"pragma no cover:{custom}")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_qphi(setting=setting)
    elif custom=="PLL":
        kwargs_likelihood["custom_logL_addition"] = logL_Prior_Lens_Light_Mltf(multifilter_setting=multifilter_setting)
    elif custom=="combined":
        #pragma=no cover
        raise RuntimeError("the logL_ellipticiy_combined has to be implemented")
        #kwargs_likelihood["custom_logL_addition"] = logL_ellipticity_combined()
    
    else:
        raise RuntimeError("The 'custom' string parameter must be :'aligned','qphi' or 'combined'(WIP). Not "+str(custom))
    return kwargs_likelihood
