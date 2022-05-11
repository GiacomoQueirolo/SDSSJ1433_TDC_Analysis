#SIMPLE COPY FROM settings_f105w_FINAL_wLL.py
#Initial data parameter settings for f105w with lens light subtracted
#####################################################################
from numpy import array, inf, nan, pi,ones_like
from conversion import conv_5 
import os
 
#Comment to insert in the backup file
comments="Change wrt to settings_f105w_FINAL_VII_ws: same lens prior for theta_E, e1,e2, theta_E_pert and shear parameters "   

R_sersic_lens = 0.11
sg_R_sersic_lens = 0.05

n_sersic_lens = 1.36 
sg_n_sersic_lens= 0.5

x_source = 0.466 
y_source = -1.74
sg_x_source = 0.02
sg_y_source = 0.02


#Modify here to indicate the image and the psf files and the model of lens light
data_path="./data/F105W/"
image_name="F105W_qso_lens_cut.fits"
err_name="e.F105W_qso_lens_cut.fits"

psf_name="big_sub_psf_0001_seA_F105W.fits"
err_psf="e.big_sub_psf_0001_seA_F105W.fits"

pssf = 1 #point source supersampling factor -> how much the PSF was supersampled
ssks = 17#it should be the size of the non-	supersampled kernel in pixel
#lens_light_model_name ="model_f814w_large.fits"
lens_light_model_name= None
#bool_bin = True #to bin the psf in the code
bool_bin=False
fixed_mag =[False] #Fixed magnification (one element per image)
sub = False #subtracted light
filter_name = "f105w"
#"cosmetics": coutout for the printed images
v_min,v_max=-1,1
e_v_min,e_v_max =-1,0
res_min, res_max= -4,4

#Coordinates
############
#CD matrix
# from degree to arcsec

CD1_1   = 1.04151273131054E-05*3600  
CD1_2   = 3.40685456789135E-05*3600                                                    
CD2_1   = 3.40685456789135E-05*3600                                                    
CD2_2   = -1.0415127313105E-05*3600  

from numpy import array as arr
transform_pix2angle = arr([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
ra_at_xy_0,dec_at_xy_0 = 0,0
    
# Parameter for the lens
center_x_pix,center_y_pix = 57.3934,57

#For the perturber
x_pert_pix,y_pert_pix = 54,41
# Parameters for the source
# order A (top), C (right), B (bottom),D (left)
x_im_pix = [70.4763, 52.8799, 42.4159, 58.7773]
y_im_pix = [45.4564, 44.7305, 54.0178, 65.5169]
z_source = 2.737
z_lens= 0.407

#Settings for cut image
edge = None

#MOD_SHIFT 
# in ra,dec considering the rotation and deformation of the reference system
mod_fr = "We consider the absolute frame of reference"
x_image, y_image, x_pert, y_pert, center_x, center_y, pix_scale ,ra_at_xy_0,dec_at_xy_0= \
	conv_5 (x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,transform_pix2angle,ra_at_xy_0,dec_at_xy_0)


#Mask coord
###########
x_mask, y_mask =[60,83,54,93,center_x_pix],[24,48,26,38,center_y_pix] # coord in pixel in the HST image 
#radius of the mask for the object(s) in pixels
r_mask = [3,3,1,2,2]
x_mask_conv,y_mask_conv = array(x_mask)-1,array(y_mask)-1

new_mask="we consider a new masking method"

#Model parameters
#################

#from f105w
theta_E = 1.645 
min_theta_E = 1.4
max_theta_E = 1.7
sig_theta_E = 0.5

# from f814 (with small variations)
e1     = 0.19
min_e1 = -.01
max_e1 = .4
sig_e1 = .1

e2     = -0.07
min_e2 = -0.3
max_e2 = 0.2
sig_e2 = .1
#Lens mass profile
fixed_lens = []
kwargs_lens_init = []
kwargs_lens_sigma = []
kwargs_lower_lens = []
kwargs_upper_lens = []

fixed_lens.append({})
kwargs_lens_init.append({'center_x': center_x, 
			'center_y':center_y,
            'theta_E': theta_E,
			'e1': e1,
		 	'e2': e2})

kwargs_lens_sigma.append({'center_x': .9*pix_scale, 
			            'center_y':.9*pix_scale,
                        'theta_E': sig_theta_E, 
                        'e1': sig_e1, 
                        'e2': sig_e2})

kwargs_lower_lens.append({'center_x': center_x-4.*pix_scale,
                          'center_y':center_y-2.*pix_scale,
                          'theta_E': min_theta_E, 
                          'e1': min_e1, 
                          'e2': min_e2})

kwargs_upper_lens.append({'center_x': center_x+2.*pix_scale, 
                          'center_y':center_y+2.*pix_scale,
                          'theta_E': max_theta_E, 
                          'e1': max_e1, 
                          'e2': max_e2})
                   
# from f814w       
# Perturber mass profile
fixed_lens.append({})
kwargs_lens_init.append({'theta_E': 0.147 ,
                         'center_x': x_pert, 
                         'center_y': y_pert})

kwargs_lens_sigma.append({'theta_E': 0.05 , 
			              'center_x': 6.*pix_scale,
                          'center_y': 6.*pix_scale})

kwargs_lower_lens.append({'theta_E': 0.,
			              'center_x': x_pert-7*pix_scale,
                          'center_y': y_pert-7*pix_scale})

kwargs_upper_lens.append({'theta_E': 0.35, 
		               	  'center_x': x_pert+7*pix_scale,
                          'center_y': y_pert+7*pix_scale})

# from f814w and f160
#Shear
fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
kwargs_lens_init.append({'gamma_ext' : 0.112, 'psi_ext': -0.41*pi})
kwargs_lens_sigma.append({'gamma_ext': 0.05, 'psi_ext': pi/50.})
kwargs_lower_lens.append({'gamma_ext': 0.00, 'psi_ext':-0.5*pi})
kwargs_upper_lens.append({'gamma_ext': 0.40, 'psi_ext':-.2*pi})
                          

#Source light parameters
fixed_source = [{}]
kwargs_source_init = [{'amp' : 0,'R_sersic': .47 ,
				 'n_sersic': 1.1, 
				 'e1': -0.1, 
				 'e2': -0.02, 
				 'center_x': x_source, 
				 'center_y': y_source}]
kwargs_source_sigma = [{'amp': 0,'R_sersic': 0.2, 
				'n_sersic': 1., 
				'e1': 0.05 , 
				'e2': 0.05,
				'center_x': sg_x_source*3, 
				'center_y': sg_y_source*3}]
kwargs_lower_source = [{'amp': 0,'R_sersic': .1, 
				'n_sersic': 0.1, 
				'e1': -.23, 
				'e2': -0.2, 
				'center_x': x_source - (sg_x_source*18), 
				'center_y': y_source - (sg_y_source*18)}]

kwargs_upper_source = [{'amp': 0,'R_sersic': 1.5, 
				'n_sersic': 6., 
				'e1': 0.0, 
				'e2': 0.05, 
				'center_x': x_source + (sg_x_source*18), 
				'center_y': y_source + (sg_y_source*18)}]

#Ps light parameters 
fixed_ps = [{}]
lum_ps = array([0,0,0,0])
kwargs_ps_init = [{'ra_image': x_image, 
		'dec_image': y_image,
		"point_amp": lum_ps}]
bound_pix_ps=2.5
kwargs_ps_sigma = [{'ra_image' : 1*pix_scale*ones_like(x_image), 
                    'dec_image': 1*pix_scale*ones_like(x_image),"point_amp":lum_ps  }]
kwargs_lower_ps = [{'ra_image' : x_image-bound_pix_ps*pix_scale, 
                    'dec_image': y_image-bound_pix_ps*pix_scale,"point_amp": lum_ps}]
kwargs_upper_ps = [{'ra_image' : x_image+bound_pix_ps*pix_scale, 
                    'dec_image': y_image+bound_pix_ps*pix_scale,"point_amp": lum_ps}]
#Lens light parameters
fixed_lens_light = []
kwargs_lens_light_init = []
kwargs_lens_light_sigma = []
kwargs_lower_lens_light = []
kwargs_upper_lens_light = []
#Lens light parameters
fixed_lens_light = [{}] 
kwargs_lens_light_init =[ {'amp': 0.,'R_sersic':2.1,
				'n_sersic':4.3, 
		 		'e1': 0.26 ,
		 		'e2': 0.05,
		 		'center_x': center_x, 
		 		'center_y': center_y}]
kwargs_lens_light_sigma=[{'amp': 0.,'R_sersic':.2, 
				'n_sersic':0.5, 
				'e1': .1,
				'e2': 0.1,
				'center_x': 0.8*pix_scale, 
				'center_y': .8*pix_scale}]
kwargs_lower_lens_light=[{'amp': 0.,'R_sersic': 1.,
				'n_sersic':1.,
			 	'e1': .15 ,
			 	'e2': -.3,
				'center_x': center_x-2.*pix_scale,
				'center_y':center_y-4.*pix_scale,}]
kwargs_upper_lens_light = [{'amp':0.,'R_sersic':5, 
				'n_sersic':6,
			 	'e1': 0.5 ,
			 	'e2': 0.1,
				'center_x': center_x+2.*pix_scale,
				'center_y':center_y+2.*pix_scale}]
# Perturber light profile
bound_pix_pert=7 #5
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					'n_sersic':n_sersic_lens , 
					'center_x':x_pert,
					'center_y':y_pert})
kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					'n_sersic': sg_n_sersic_lens*3, 
					'center_x':6.*pix_scale,
					'center_y':6.*pix_scale}) 
kwargs_lower_lens_light.append({'amp':0.,'R_sersic':.1, 
					'n_sersic': 1,
	                                'center_x':x_pert - bound_pix_pert*pix_scale,
	                                'center_y':y_pert - bound_pix_pert*pix_scale}) 
kwargs_upper_lens_light.append({'amp':0.,'R_sersic':3, 
					'n_sersic': 8,
	                                'center_x':x_pert + bound_pix_pert*pix_scale,
	                                'center_y':y_pert + bound_pix_pert*pix_scale}) 

#Uniform light profile
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp':1.5})
kwargs_lens_light_sigma.append({'amp':3})
kwargs_lower_lens_light.append({'amp':0.})
kwargs_upper_lens_light.append({'amp':10.}) 

#MOD_PROD
def produce_kwargs_result(samples_mcmc,param_mcmc,i):
	from produce_kwargs_result import with_main_lens_light
	return with_main_lens_light(samples_mcmc,param_mcmc,i)
	
def str_param(param):
    from get_param_title import standard
    return standard(param)
