from Utils.tools import *

class regdiff_set():
    def __init__(self,name,lcs=None,pd=2.,covkernel="matern",pow=1.5,\
                 errscale=1,amp=1,scale=200,\
                 verbose=True, method="weights"):
        """
        :param method: optimisation method. Choose between "weights" and "simple" (default : "weights")
        :type method: str
        :param pd: the point density, in points per days.
        :type pd: float
        :param covkernel: Choose between "matern","RatQuad" and "RBF". See scikit GP documentation for details
        :type covkernel: str
        :param pow: float, exponent coefficient of the covariance function
        :type pow: float
        :param amp: initial amplitude coefficient of the covariance function
        :type amp: float
        :param scale: float, initial characteristic time scale
        :type scale: amp
        :param errscale: additional scaling of the photometric error
        :type errscale: float

        """
        self.set_name  = str(name)
        self.lcs       = lcs
        self.pd        = pd
        self.covkernel = covkernel
        self.pow       = pow
        self.errscale  = errscale
        self.amp       = amp
        self.scale     = scale
        self.verbose   = verbose
        self.method    = method
        #self.rgd_set_dir = self.def_dir(name=name)        
    def get_kw(self):
        return self.__dict__
    def __str__(self):
        _str = r""
        kw   = self.get_kw()
        for k in kw:
            _str+=str(k)+" "+str(kw[k])+","
    

from Utils.tools import *
from TD_analysis.pycs_get_res import *
from TD_analysis.stnd_handling_data import get_ref_index


def find_sim(config):
    ref_config = get_config(config.refer_config)
    group_res = get_group_res(ref_config)
    best_res = group_res[get_ref_index(group_res)]
    sim_path_best = str(best_res.sim_path)
    mocks_best = sim_path_best+"/sims_mocks_n800t10_PS/"

    sim_dir    = config.simu_directory
    os.symlink(sim_path_best+)