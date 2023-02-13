from astropy import units as u
from Utils.tools import print_res_w_err
import numpy as np

class H0_Res():
    def __init__(self,H0,sigma):
        self.unit = u.km/(u.second*u.Mpc)
        self.H0 = H0
        if type(sigma) is list:
            self.sigma_min = sigma[0]
            self.sigma_max = sigma[1]
            self.sigma     = np.mean(sigma) 
        else:
            self.sigma_min = sigma
            self.sigma_max = sigma
            self.sigma     = sigma 
    def __str__(self):
        return print_res_w_err(self.H0,self.sigma)+" "+str(self.unit)
    
# planck  = 67.4 +- 0.5 km/s/Mpc
h0planck = H0_Res(H0=67.4,sigma=0.5)
# H0LiCoW = 73.3 -1.8+1.7  km/s/Mpc
h0licow  = H0_Res(H0=73.3,sigma=[1.8,1.7])

