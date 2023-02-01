import numpy as np
from conversion import e12_from_qphi as e12
class lens_prior(object):
    #main lens
    # changed e1,e2 defined from q,phi
    theta_E     = 1.65 
    min_theta_E = 1.5
    max_theta_E = 1.8
    sig_theta_E = 0.5

    # mashup: e1 from same_prior_ellipt_I
    # e2 from  same_prior_ellipt_II 
    # e2 -> extended even more, hoping that the logL change will be enough
    # same for e1, more generic prior
    ################################################
    e1             = 0.15  #0.1322
    min_e1         = 0.0   #0.0192
    max_e1         = 0.3   #0.2716
    sig_e1         = 0.2   #0.05
    
    e2             = -0.00 #0422#-0.0411
    min_e2         = -0.12 #1807#-0.1033
    max_e2         =  0.12 #0.0090#0306#-0.0053
    sig_e2         =  0.20 #0.05  # 0.02
    ################################################
    
    #pert
    theta_pert     = 0.15
    min_theta_pert = 0.00
    max_theta_pert = 0.35
    sig_theta_pert = 0.05
    
    #shear
    gamma_ext     = 0.11
    min_gamma_ext = 0.00
    max_gamma_ext = 0.40
    sig_gamma_ext = 0.05
    
    psi_ext     = -0.41*np.pi
    min_psi_ext = -0.5*np.pi 
    max_psi_ext = -0.2*np.pi
    sig_psi_ext = np.pi/50.

class lens_PEMD_prior(lens_prior):
    # free gamma
    gamma     = 2.0
    sig_gamma = 0.3
    min_gamma = 1.5
    max_gamma = 2.5 
