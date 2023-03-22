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

    def _get_other_H0(self,other):
        if isinstance(other,H0_Res):
            H0_prime=other.H0
        elif isinstance(other,float) or isinstance(other,int):
            H0_prime=other 
        else:
            raise TypeError(f"{other} not recognised, type must be H0_Res, float or int, not {type(other)}")
        return H0_prime

    def __eq__(self, other):
        H0_prime = self._get_other_H0(other)
        return self.H0==H0_prime
    
    def __ge__(self, other):
        H0_prime = self._get_other_H0(other)
        return self.H0>=H0_prime
    
    def __gt__(self, other):
        H0_prime = self._get_other_H0(other)
        return self.H0>H0_prime
    
    def __lt__(self, other):
        H0_prime = self._get_other_H0(other)
        return self.H0<H0_prime
    
    def __le__(self, other):
        H0_prime = self._get_other_H0(other)
        return self.H0<=H0_prime
    
    def __float__(self):
        return self.H0

    def __add__(self,other):
        H0_prime = self._get_other_H0(other)
        return self.H0+H0_prime
        
    def __sub__(self,other):
        H0_prime = self._get_other_H0(other)
        return self.H0-H0_prime



# planck  = 67.4 +- 0.5 km/s/Mpc
h0planck = H0_Res(H0=67.4,sigma=0.5)
# H0LiCoW = 73.3 -1.8+1.7  km/s/Mpc
h0licow  = H0_Res(H0=73.3,sigma=[1.8,1.7])

