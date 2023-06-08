import numpy as np
# copy for same_prior_redone.py but with better structure
# get_** must be independent of the values, which must be given as input

class lens_prior(object):
    z_source = 2.737
    z_lens   = 0.407

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

    bound_pix_mainlens = 4.
    bound_pix_pert     = 4.


class lens_PEMD_prior(lens_prior):
    # free gamma
    gamma     = 2.0
    sig_gamma = 0.3
    min_gamma = 1.5
    max_gamma = 2.5 


def get_lens_params(center_x, center_y, x_pert, y_pert, pix_scale, CP):
    prior_lens        = lens_prior()
    bound_mainlens    = prior_lens.bound_pix_mainlens*pix_scale
    sigma_mainlens    = .9

    bound_pert        = prior_lens.bound_pix_pert*pix_scale
    sigma_pert        = 2.
    
    fixed_lens        = []
    kwargs_lens_init  = []
    kwargs_lens_sigma = []
    kwargs_lower_lens = []
    kwargs_upper_lens = []

    fixed_lens.append({})
    kwargs_lens_init.append({'theta_E':  prior_lens.theta_E,
                             'center_x': center_x, 
		                     'center_y': center_y,
		                     'e1':       prior_lens.e1,
	         	             'e2':       prior_lens.e2})

    kwargs_lens_sigma.append({'theta_E':  prior_lens.sig_theta_E, 
                              'center_x': sigma_mainlens*pix_scale, 
		                      'center_y': sigma_mainlens*pix_scale,
                              'e1':       prior_lens.sig_e1, 
                              'e2':       prior_lens.sig_e2})

    kwargs_lower_lens.append({'theta_E':  prior_lens.min_theta_E, 
                              'center_x': center_x-bound_mainlens,
                              'center_y': center_y-bound_mainlens,
                              'e1':       prior_lens.min_e1, 
                              'e2':       prior_lens.min_e2})

    kwargs_upper_lens.append({'theta_E':  prior_lens.max_theta_E, 
		                      'center_x': center_x+bound_mainlens, 
                              'center_y': center_y+bound_mainlens,
                              'e1':       prior_lens.max_e1, 
                              'e2':       prior_lens.max_e2})
    if CP:
        prior_lens=lens_PEMD_prior()
        kwargs_lens_init[0]["gamma"]  = prior_lens.gamma
        kwargs_lens_sigma[0]["gamma"] = prior_lens.sig_gamma
        kwargs_lower_lens[0]["gamma"] = prior_lens.min_gamma
        kwargs_upper_lens[0]["gamma"] = prior_lens.max_gamma

    # Perturber mass profile


    fixed_lens.append({})
    kwargs_lens_init.append({'theta_E' : prior_lens.theta_pert,
                             'center_x': x_pert, 
                             'center_y': y_pert})

    kwargs_lens_sigma.append({'theta_E' : prior_lens.sig_theta_pert, 
    			              'center_x': sigma_pert*pix_scale,
                              'center_y': sigma_pert*pix_scale})

    kwargs_lower_lens.append({'theta_E' : prior_lens.min_theta_pert,
		                      'center_x': x_pert-bound_pert,
                              'center_y': y_pert-bound_pert})

    kwargs_upper_lens.append({'theta_E' : prior_lens.max_theta_pert,
		                      'center_x': x_pert+bound_pert,
                              'center_y': y_pert+bound_pert})
                              
    #Shear
    fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
    kwargs_lens_init.append({'gamma_ext' : prior_lens.gamma_ext,     'psi_ext': prior_lens.psi_ext})
    kwargs_lens_sigma.append({'gamma_ext': prior_lens.sig_gamma_ext, 'psi_ext': prior_lens.sig_psi_ext})
    kwargs_lower_lens.append({'gamma_ext': prior_lens.min_gamma_ext, 'psi_ext': prior_lens.min_psi_ext})
    kwargs_upper_lens.append({'gamma_ext': prior_lens.max_gamma_ext, 'psi_ext': prior_lens.max_psi_ext})
    
    lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
    
    return lens_params
    
def get_ps_params(x_image,y_image,pix_scale,fixed_pos=False):
    #Ps light parameters

    lum_ps   = np.array([0,0,0,0])
    bound_pix_ps   = 2.# 1.8 #.6
    bound_ps       = bound_pix_ps*pix_scale
    sigma_im       = 1.
    if not fixed_pos:
        fixed_ps = [{}]    
        kwargs_ps_init  = [{'ra_image' : x_image, 'dec_image':y_image,"point_amp": lum_ps}]
        kwargs_ps_sigma = [{'ra_image' : sigma_im*pix_scale*np.ones_like(x_image), 
                            'dec_image': sigma_im*pix_scale*np.ones_like(x_image),"point_amp": lum_ps}]
        kwargs_lower_ps = [{'ra_image' : x_image-bound_ps, 'dec_image': y_image-bound_ps,"point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : x_image+bound_ps, 'dec_image': y_image+bound_ps,"point_amp": lum_ps}]
    else:
        fixed_ps = [{'ra_image' : x_image, 'dec_image':y_image}]    
        kwargs_ps_init  = [{"point_amp": lum_ps}]
        kwargs_ps_sigma = [{"point_amp": lum_ps}]
        kwargs_lower_ps = [{"point_amp": lum_ps}]
        kwargs_upper_ps = [{"point_amp": lum_ps}]
    ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
    return ps_params
    
def get_lens_light_params(center_x,center_y,x_pert,y_pert,pix_scale,SUB):
    prior_lens        = lens_prior()
    bound_mainlens    = prior_lens.bound_pix_mainlens*pix_scale
    bound_pert        = prior_lens.bound_pix_pert*pix_scale
    #Lens light parameters
    fixed_lens_light        = []
    kwargs_lens_light_init  = []
    kwargs_lens_light_sigma = []
    kwargs_lower_lens_light = []
    kwargs_upper_lens_light = []
    if not SUB:
        #Lens light parameters
        fixed_lens_light = [{}] 
        kwargs_lens_light_init =[ {'amp': 0.,
                        'R_sersic':1.5,
				        'n_sersic':4.5, 
		         		'e1': 0.26 , #already optimised
		         		'e2': 0.05 , #already optimised
		         		'center_x': center_x, 
		         		'center_y': center_y}]
        kwargs_lens_light_sigma=[{'amp': 0.,'R_sersic':.2, 
				        'n_sersic':0.2,  
				        'e1': .1,
				        'e2': 0.1,
				        'center_x': 0.8*pix_scale, 
				        'center_y': .8*pix_scale}]
        kwargs_lower_lens_light=[{'amp': 0.,'R_sersic': 1.,
				        'n_sersic':2.5,
			         	'e1': .15 ,
			         	'e2': -.3,
				        'center_x': center_x-bound_mainlens,
				        'center_y': center_y-bound_mainlens}]
        kwargs_upper_lens_light = [{'amp':0.,'R_sersic':5, 
				        'n_sersic':6,
			         	'e1': 0.5 ,
			         	'e2': 0.1,
				        'center_x': center_x+bound_mainlens,
				        'center_y': center_y+bound_mainlens}]
    # Perturber light profile
    fixed_lens_light.append({})
    kwargs_lens_light_init.append({'amp' :0.,'R_sersic': 0.363 ,
				                'n_sersic': 4., 
				                'center_x':x_pert,
				                'center_y':y_pert})
    kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': .3,  
				                    'n_sersic': 1.5, 
				                    'center_x': 2.*pix_scale,
				                    'center_y': 2.*pix_scale}) 
    kwargs_lower_lens_light.append({'amp':0.,'R_sersic': .1, 
				                    'n_sersic': 1.,
                                    'center_x':x_pert - bound_pert,
                                    'center_y':y_pert - bound_pert}) 
    kwargs_upper_lens_light.append({'amp':0.,'R_sersic': 2., 
				                     'n_sersic': 5.,
                                     'center_x':x_pert + bound_pert,
                                     'center_y':y_pert + bound_pert}) 

    #Uniform light profile
    fixed_lens_light.append({})
    kwargs_lens_light_init.append( {'amp':1.5})
    kwargs_lens_light_sigma.append({'amp':3})
    kwargs_lower_lens_light.append({'amp':0.})
    kwargs_upper_lens_light.append({'amp':10.}) 
    
    lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]
    
    return lens_light_params
    
    
    
def get_source_params(pixscale):
    x_source     = 0.5
    y_source     = -1.7
    sig_cent     = 0.05
    
    fixed_source       = [{}]
    kwargs_source_init = [{'amp' : 0,'R_sersic': 0.3 ,
			                'n_sersic': 1.4, 
			                'center_x': x_source, 
			                'center_y': y_source}]
    kwargs_source_sigma = [{'amp': 0,'R_sersic': 0.06, 
			                'n_sersic': 1, 
			                'center_x': sig_cent, 
			                'center_y': sig_cent}]
    kwargs_lower_source = [{'amp': 0,'R_sersic': 0.1, 
			                'n_sersic': 1., 
			                'center_x': x_source - (sig_cent), 
			                'center_y': y_source - (sig_cent)}]
    kwargs_upper_source = [{'amp': 0,'R_sersic': 0.5, 
			                'n_sersic': 7., 
			                'center_x': x_source + (sig_cent), 
			                'center_y': y_source + (sig_cent)}]
    source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    return source_params
