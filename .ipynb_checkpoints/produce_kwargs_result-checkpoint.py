#MOD_PROD
# list of possible produce_kwargs_result functions
# NOTE: kwargs_result is actually useless for the produce_kwargs_result function itself - rework them all
import numpy as np

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

    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
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
    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
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

    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
    return kwargs_result_i

def newProfile_with_main_lens_light(samples_mcmc,param_mcmc,i):
    """
    for the model with PEMD profile with main lens light modelled (for F160W)
    
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

    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
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

    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
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

    kwargs_result_i["kwargs_ps"] = [{"ra_image":np.array([samples_mcmc[i][param_mcmc.index('ra_image'):param_mcmc.index('ra_image')+4]]),
    'dec_image':np.array([samples_mcmc[i][param_mcmc.index('ra_image')+4:param_mcmc.index('ra_image')+8]]) }]
    
    return kwargs_result_i
