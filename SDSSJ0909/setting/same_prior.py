from numpy import pi

class lens_prior(object):
    #main lens
    theta_E     = 5.2
    min_theta_E = 4.
    max_theta_E = 7.
    sig_theta_E = 0.5

    e1     =  0.
    min_e1 = -0.7
    max_e1 =  0.7
    sig_e1 =  0.1

    e2     =  0.0
    min_e2 = -0.7
    max_e2 =  0.7
    sig_e2 =  0.1
    
    #shear
    gamma_ext     = 0.1
    min_gamma_ext = 0.00
    max_gamma_ext = 0.50
    sig_gamma_ext = 0.1
    
    psi_ext     = 0
    min_psi_ext = -1.*pi 
    max_psi_ext = +1*pi
    sig_psi_ext = pi/10.

class lens_PEMD_prior(lens_prior):
    # free gamma
    gamma     = 2.0
    sig_gamma = 0.3
    min_gamma = 1.5
    max_gamma = 2.5 
