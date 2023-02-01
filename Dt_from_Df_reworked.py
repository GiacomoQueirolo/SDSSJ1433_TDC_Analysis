#!/usr/bin/env python
# coding: utf-8

# # From Dt to Df and viceversa

# In[1]:


import copy
import pickle
import numpy as np
from astropy import units as u
import scipy.integrate as integrate
from astropy import constants as const
from astropy.cosmology import default_cosmology
   
def get_cosmo_prm(setting):
    z_source  = setting.z_source
    z_lens    = setting.z_lens

    
    cosmo_dd  = cosmo.angular_diameter_distance(z_lens).to("kpc")   #kpc
    cosmo_ds  = cosmo.angular_diameter_distance(z_source).to("kpc") #kpc
    cosmo_dds = cosmo.angular_diameter_distance_z1z2(z1=z_lens,z2=z_source).to("kpc") #kpc
    


# In[30]:


def k_analytical(setting,cosmo=None):
    if not cosmo:
        cosmo  = default_cosmology.get()
    z_d = setting.z_lens
    z_s = setting.z_source
    Dd  = cosmo.angular_diameter_distance(z_d)
    Ds  = cosmo.angular_diameter_distance(z_s)
    Dds = cosmo.angular_diameter_distance_z1z2(z1=z_d,z2=z_s)
    H0  = cosmo.H0
    k_test = (1+z_d)*Dd*Ds*H0/(const.c*Dds)
    return k_test.to("").value

def H_0(D_phi,D_t,k):
    D_t   *=u.day.to("s")             # convert from day in sec
    D_phi *=u.arcsec.to("degree")**2  # convert from arcsec^2 to degree^2
    D_phi *=u.degree.to("rad")**2     # convert from degree^2 to rad^2 (non-dim)
    h_0    = k*D_phi/D_t              # in unit of 1/sec (~frequency)
    H_0    = h_0*Mpc/km               # convert to km/(s*Mpc)
    return H_0

def sig_H0(D_phi,D_t,k,sig_phi=0,sig_Dt=0):
    D_t     *=u.day.to("s")             # convert from day in sec
    sig_Dt  *=u.day.to("s")             # convert from day in sec
    D_phi   *=u.arcsec.to("degree")**2  # convert from arcsec^2 to degree^2
    D_phi   *=u.degree.to("rad")**2     # convert from degree^2 to rad^2 
    sig_phi *=u.arcsec.to("degree")**2  # convert from arcsec^2 to degree^2
    sig_phi *=u.degree.to("rad")**2     # convert from degree^2 to rad^2 
    #sig_ho = k(z_l,z_s) * np.sqrt((sig_phi/D_t)**2 + (D_phi*sig_Dt/(D_t**2))**2) #unit of 1/sec
    sig_ho  = (k/(D_t**2)) * np.sqrt((D_t*sig_phi)**2 + (D_phi*sig_Dt)**2) #unit of 1/sec (simplification)
    sig_ho  = sig_ho*Mpc/km             # unit of km/(s*Mpc)
    return sig_ho


def Dt_XY (Df_XY,H0,k): 
    """
    # Gives the time delay [days] of Y wrt X, obtained from Delta fermat, which is
    # the difference in fermat potential at position Y wrt position X (usually X=A)
    # Param:
    ########
    # Df_XY: float or array
    #     differences in fermat pot btw position Y and X in arcsec^2
    # H0 : float
    #    H0 in km/s/Mpc
    # k  : float
    #     = (1+z_lens)*H0*Dd*Ds/(c*Dds) (non-dimensional, small cosmological denpendence)
    # Return:
    #########
    # Dt_XY: float or array
    #     time delay btw position Y and X in days
    """
    Df_XY_  = copy.deepcopy(Df_XY)
    H0_     = copy.deepcopy(H0)
    Df_XY_ *= u.arcsec.to("rad")**2       # convert from arcsec^2 to rad^2
    H0_    *= u.km.to("Mpc")              # convert from km/s/Mpc to 1/s
    Dt_XY   = k*Df_XY_/H0_
    Dt_XY  *= u.second.to("d")            # convert from sec to day
    return Dt_XY

def Df_XY (Dt_XY,H0,k): 
    """
    # Gives the Delta fermat [arcsec^2] of Y wrt X, obtained from Delta time, which is
    # the difference in arrival time at position Y wrt position X (usually X=A)
    # Param:
    ########
    # Dt_XY: float or array
    #     differences in arrival time btw position Y and X in days
    # H0 : float
    #    H0 in km/s/Mpc
    # k  : float
    #     = (1+z_lens)*H0*Dd*Ds/(c*Dds) (non-dimensional, small cosmological denpendence)
    # Return:
    #########
    # Df_XY: float or array
    #     differences in fermat pot btw position Y and X in arcsec^2
    """
    Dt_XY_  = copy.deepcopy(Dt_XY)
    H0_     = copy.deepcopy(H0)
    Dt_XY_ *= u.day.to("s")            # convert from day to sec
    H0_    *= u.km.to("Mpc")           # convert from km/s/Mpc to 1/s
    Df_XY   = H0_*Dt_XY_/k
    Df_XY  *= u.rad.to("arcsec")**2    # convert from rad^2 to arcsec^2
    return Df_XY

def test_DtDf_XY(n=1000):
    class sett():
        def __init__(self):
            self.z_source = 2.7
            self.z_lens   = 0.4
    k  = k_analytical(sett())
    for Dt in np.random.uniform(1,100,int(n/10)):
        for h0 in  np.random.uniform(50,100,int(n/10)):
            Dt_ = copy.copy(Dt)
            for i in range(n):
                Df_ = Df_XY(Dt_,h0,k)
                Dt_ = Dt_XY(Df_,h0,k)
            if abs(Dt-Dt_)>.1:
                err_message = "the inversion of Dt in Df is diverging for H0=",h0,"km/s/Mpc and Dt=",Dt,"d"
                raise RuntimeError(err_message)
    return 0

def cov_Df(cov_Dt,H0,k): 
    """
    # Gives the covariance of  Delta fermat [(arcsec^2)^2] of Y wrt X, 
    # obtained from cov. Delta time, which is the cov of the difference
    # in arrival time at position Y wrt position X (usually X=A)
    # follow from https://statproofbook.github.io/P/mvn-ltt.html
    # Param:
    #######
    # cov_Dt: float or array
    #     cov. of differences in arrival time btw position Y and X in days
    # H0 : float
    #    H0 in km/s/Mpc
    # k  : float
    #     = (1+z_lens)*H0*Dd*Ds/(c*Dds) (non-dimensional, small cosmological denpendence)
    # Return:
    #########
    # cov_Df: float or array
    #   cov. of differences in fermat pot btw position Y and X in arcsec^2
    """
    cov_Dt_  = copy.deepcopy(cov_Dt)
    H0_      = copy.deepcopy(H0)
    cov_Dt_ *= u.day.to("s")**2         # convert from day^2 to sec^2
    H0_     *= u.km.to("Mpc")           # convert from km/s/Mpc to 1/s
    cov_Df   = cov_Dt_*((H0_/k)**2)
    cov_Df  *= u.rad.to("arcsec")**4    # convert from rad^2 to arcsec^2
    return cov_Df


# In[33]:


#test_DtDf_XY(100)=Success


# In[ ]:


#########################
# old stuff
"""
#Constants
c   = const.c.to("m/s") #m/s
pc  = u.parsec.to("m")# m
Mpc = u.Mpc.to("m") #m
km  = u.km.to("m") #m
# Cosmologic paramaters 
cosmo   = default_cosmology.get()
Omega_m = cosmo.Om0   # omega matter = 0.3075
Omega_k = cosmo.Ok0   # omega curvature = 0.0
Omega_L = 1-Omega_m   # omega lambda 

# SDSSJ1433

z_s = 2.737
z_l = 0.407

def E(z,Omega_m0=Omega_m,Omega_L0=Omega_L):
    return np.sqrt((Omega_m0*((1+z)**3) + Omega_L0))  #non-dimensional

def k(z_l,z_s,Omega_m0=Omega_m,Omega_L0=Omega_L): 
    zs    = integrate.quad(lambda x: (1./E(x)) ,0,z_s)[0]
    zl    = integrate.quad(lambda x: (1./E(x)) ,0,z_l)[0]
    zl_zs = integrate.quad(lambda x: (1./E(x)) ,z_l,z_s)[0]
    return zs*zl/zl_zs #non-dimensional

# SDSSJ1433

z_s = 2.737
z_l = 0.407

def E(z,Omega_m0=Omega_m,Omega_L0=Omega_L):
    return np.sqrt((Omega_m0*((1+z)**3) + Omega_L0))  #non-dimensional

def k(z_l,z_s,Omega_m0=Omega_m,Omega_L0=Omega_L): 
    zs    = integrate.quad(lambda x: (1./E(x)) ,0,z_s)[0]
    zl    = integrate.quad(lambda x: (1./E(x)) ,0,z_l)[0]
    zl_zs = integrate.quad(lambda x: (1./E(x)) ,z_l,z_s)[0]
    return zs*zl/zl_zs #non-dimensional


def E(z,Omega_m0=Omega_m,Omega_L0=Omega_L):
    #non-dimensional
    return np.sqrt(Omega_m0*((1+z)**3) + Omega_L0)

def k(z_l,z_s,Omega_m0=Omega_m,Omega_L0=Omega_L): 
    #also = (1+z_d) * D_d D_d /D_ds
    #non-dimensional
    zs   = integrate.quad(lambda x: (1./E(z,Omega_m0=Omega_m,Omega_L0=Omega_L)) ,0,z_s)[0]
    zl   = integrate.quad(lambda x: (1./E(z,Omega_m0=Omega_m,Omega_L0=Omega_L)) ,0,z_l)[0]
    zl_zs= integrate.quad(lambda x: (1./E(z,Omega_m0=Omega_m,Omega_L0=Omega_L)) ,z_l,z_s)[0]
    return zs*zl/zl_zs 

def H_0(D_phi,D_t,z_l=z_l,z_s=z_s):
    D_t   *=24.*60.*60.               # convert from day in sec
    D_phi /=(3600.*3600.)             # convert from arcsec^2 to degree^2
    D_phi *=(np.pi*np.pi/(180.*180.)) # convert from degree^2 to rad^2 
    h_0 = k(z_l,z_s)*D_phi/D_t        # in unit of 1/sec
    H_0 = h_0*Mpc/km                  # convert to km/(s*Mpc)
    return H_0
"""

