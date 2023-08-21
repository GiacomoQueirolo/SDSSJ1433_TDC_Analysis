
# for HST_Main.py
# recreate initial mcmc step from previous PSO backup results
# so that we don't have to repeat it
import numpy as np
from lenstronomy.Util import sampling_util

from Utils.tools import *
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.get_res import load_whatever
from Data.Multiband_Data.Param_mltf import get_Param_mltf,get_sigma_mltf


@check_mltf_setting
def create_mcmc_init_mltf(multifilter_setting):
    savemcmc_path =  multifilter_setting.savemcmc_path
    try:
        try:
            pso_file  = multifilter_setting.get_savejson_name(filename="pso")
            pso_init  = load_whatever(pso_file)
            # struct: likl, pos, distances_from_best
            lkl,pos,_ = pso_init
        except FileNotFoundError:
            pso_file  = savemcmc_path+"/psobackup.json"
            pso_init  = load_whatever(pso_file)
            # struct : Name, [lkl,pos,vel], param
            _,steps,_ = pso_init
            lkl,pos,_ = steps

        # we recreate the "mc_init_sample" such at it would be at the end of the PSO run
        # in lenstronomy
        best_pso    = pos[np.array(lkl).argmax()]
        
        Prm         = get_Param_mltf(multifilter_setting)
        prm_lmt     = Prm.param_limits()
        lower_limit = prm_lmt[0]
        upper_limit = prm_lmt[1]
        
        sigma_init  = get_sigma_mltf(multifilter_setting)
        
        mc_init_sample =  sampling_util.sample_ball_truncated(best_pso,sigma_init,lower_limit,upper_limit)

        return np.array(mc_init_sample)
    except FileNotFoundError:
        print("No PSO files found, returning None mcmc initial sample")
        return None