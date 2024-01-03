
# for HST_Main.py
# recreate initial mcmc step from previous PSO backup results
# so that we don't have to repeat it
import numpy as np
from lenstronomy.Util import sampling_util

from Utils.tools import *
from Utils.get_res import load_whatever
from Data.Param import get_Param,get_sigma


@check_setting
def create_mcmc_init(setting,walkerRatio=10,backup_path="./backup_results/"):
    savemcmc_path = get_savemcmcpath(setting,backup_path=backup_path)
    try:
        try:
            pso_file  = save_json_name(setting,savemcmc_path,filename="pso")
            pso_init  = load_whatever(pso_file)
            # struct: likl, pos, distances_from_best
            lkl,pos,_ = pso_init
        except FileNotFoundError:
            pso_file  = savemcmc_path+"/psobackup.json"
            pso_init  = load_whatever(pso_file)
            #old:# struct : Name, [lkl,pos,vel], param
            # _,steps,_ = pso_init
            # lkl,pos,_ = steps
            lkl,pos,_ = pso_init
            
        # we recreate the "mc_init_sample" such at it would be at the end of the PSO run
        # in lenstronomy
        best_pso    = pos[np.array(lkl).argmax()]
        
        Prm         = get_Param(setting)
        prm_lmt     = Prm.param_limits()
        lower_limit = prm_lmt[0]
        upper_limit = prm_lmt[1]
        
        sigma_init  = Prm.kwargs2args(**get_sigma(setting))
        
        num_param   = Prm.num_param()[0]
        n_walkers   = num_param * walkerRatio        
        mc_init_sample =  sampling_util.sample_ball_truncated(best_pso,sigma_init,lower_limit,upper_limit,size=n_walkers)

        return np.array(mc_init_sample)
    except FileNotFoundError:
        print("No PSO files found, returning None mcmc initial sample")
        return None
