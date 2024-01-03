#MOD_PROD
# list of possible produce_kwargs_result functions
# NOTE: kwargs_result is actually useless for the produce_kwargs_result function itself - rework them all
import numpy as np


def _get_kwargs_ps(samples_mcmc,param_mcmc,i):
    return [{"ra_image":np.array(samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]),
    'dec_image':np.array(samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]) }]
    
def standard(samples_mcmc,param_mcmc,i):
    """
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    
    return kwargs_result_i


def with_main_lens_light(samples_mcmc,param_mcmc,i):
    """
    for the model where the main lens light profile is not subtracted 
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light1')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light1')]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]
    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    return kwargs_result_i

def stnd_FL(samples_mcmc,param_mcmc,fixed_pos,i):
    """
    standard model BUT with fixed lenseS positions
    
    input the sample, which parameter to take and which i-step of the sample, along with the dict fixed_pos:
    fixed_pos = {"center_x_lens0":,"center_y_lens0":,"center_x_lens1":,"center_y_lens1":}
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":fixed_pos['center_x_lens0'],
                    "center_y":fixed_pos['center_y_lens0']},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":fixed_pos['center_x_lens1'],
                    "center_y":fixed_pos['center_y_lens1']
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]
    kwargs_result_i["kwargs_ps"] = _get_kwargs_ps(samples_mcmc,param_mcmc,i)    
    return kwargs_result_i
    
def stnd_FL_wo_pert(samples_mcmc,param_mcmc,fixed_pos,i,SUB):
    """
    standard model BUT with fixed MAIN lense positions -> ONLY MAIN LENS
    
    input the sample, which parameter to take and which i-step of the sample, along with the dict fixed_pos:
    fixed_pos = {"center_x_lens0":,"center_y_lens0":}
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement
    """
    if SUB:
        # without Main Lens Light:
        # the pert. pos. is coupled with its lens light pos
        # and it's then the 1st lens light profile (python start from 0):
        lens_light = "lens_light0" 
    else:
        # with Main Lens Light:
        # the pert. pos. is coupled with its lens light pos
        # and it's then the 2nd lens light profile (python start from 0):
        lens_light = "lens_light1" 
        
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":fixed_pos['center_x_lens0'],
                    "center_y":fixed_pos['center_y_lens0']},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index("center_x_"+lens_light)],
                    "center_y":samples_mcmc[i][param_mcmc.index("center_y_"+lens_light)]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]
    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    return kwargs_result_i  
    
def stnd_FS(samples_mcmc,param_mcmc,fixed_val,i,SUB):
    """
    standard model BUT with fixed SHEAR
    
    input the sample, which parameter to take and which i-step of the sample, along with the dict fixed_pos:
    fixed_val = {'ra_0': ra_0, 'dec_0':  dec_0,'gamma_ext' :  gamma_ext,'psi_ext': psi_ext})
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement
    """
    if SUB:
        # without Main Lens Light:
        # the pert. pos. is coupled with its lens light pos
        # and it's then the 1st lens light profile (python start from 0):
        lens_main = "lens0" 
        lens_pert = "lens_light0" 
    else:
        # with Main Lens Light:
        # the pert. pos. is coupled with its lens light pos
        # and it's then the 2nd lens light profile (python start from 0):
        lens_main = "lens_light0" 
        lens_pert = "lens_light1" 
        
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_'+lens_main)],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_'+lens_main)]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index("center_x_"+lens_pert)],
                    "center_y":samples_mcmc[i][param_mcmc.index("center_y_"+lens_pert)]
                    },
                    {"gamma_ext":fixed_val['gamma_ext'],
                    "psi_ext":fixed_val['psi_ext']}
                    ]
    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)  
    return kwargs_result_i  

    
def newProfile(samples_mcmc,param_mcmc,i):
    """
    for the model with PEMD profile
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "gamma":samples_mcmc[i][param_mcmc.index("gamma_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    
    return kwargs_result_i

def newProfile_with_main_lens_light(samples_mcmc,param_mcmc,i):
    """
    for the model with PEMD profile with main lens light not subtracted (for F160W)
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "gamma":samples_mcmc[i][param_mcmc.index("gamma_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light1')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light1')]},
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    
    return kwargs_result_i
    
def newProfile_with_fixed_gamma(samples_mcmc,param_mcmc,gamma_fixed,i):
    """
    for the model with PEMD profile
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "gamma":gamma_fixed,
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"] = _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    
    return kwargs_result_i
    
def newProfile_noPert(samples_mcmc,param_mcmc,i):
    """
    for the model with PEMD profile where perturber light is already subtracted
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "gamma":samples_mcmc[i][param_mcmc.index("gamma_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens1')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens1')]
                    },
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"] =  _get_kwargs_ps(samples_mcmc,param_mcmc,i)
    
    return kwargs_result_i
    

def newProfile_with_MLL_and_fixed_QSO(sett,samples_mcmc,param_mcmc,i):
    """
    for the model with PEMD profile with main lens light not subtracted and qso pos fixed (for F140W testprof)
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "gamma":samples_mcmc[i][param_mcmc.index("gamma_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light1')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light1')]},
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"]                 = sett.ps_params[0]
    kwargs_result_i["kwargs_ps"][0]["ra_image"]  = sett.ps_params[2][0]["ra_image"]
    kwargs_result_i["kwargs_ps"][0]["dec_image"] = sett.ps_params[2][0]["dec_image"]
    
    return kwargs_result_i
    
def SIEProfile_with_MLL_and_fixed_QSO(sett,samples_mcmc,param_mcmc,i):
    """
    for the model with SIE profile with main lens light not subtracted and qso pos fixed (for F140W testprof)
    
    return the kwargs_result needed for this model specifically
    needed for the fermat potential measurement 
    """
    kwargs_result_i= {}
    kwargs_result_i["kwargs_lens"]= [{"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens0")],
                    "e1":samples_mcmc[i][param_mcmc.index("e1_lens0")],
                    "e2":samples_mcmc[i][param_mcmc.index("e2_lens0")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light0')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light0')]},
                    {"theta_E":samples_mcmc[i][param_mcmc.index("theta_E_lens1")],
                    "center_x":samples_mcmc[i][param_mcmc.index('center_x_lens_light1')],
                    "center_y":samples_mcmc[i][param_mcmc.index('center_y_lens_light1')]},
                    {"gamma_ext":samples_mcmc[i][param_mcmc.index('gamma_ext_lens2')],
                    "psi_ext":samples_mcmc[i][param_mcmc.index('psi_ext_lens2')]}
                    ]

    kwargs_result_i["kwargs_ps"]                 = sett.ps_params[0]
    kwargs_result_i["kwargs_ps"][0]["ra_image"]  = sett.ps_params[2][0]["ra_image"]
    kwargs_result_i["kwargs_ps"][0]["dec_image"] = sett.ps_params[2][0]["dec_image"]
    
    return kwargs_result_i
