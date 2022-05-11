#Initial data parameter settings for f140w with lens light subtracted
#####################################################################
from numpy import array, inf, nan, pi,ones_like
from conversion import conv_5
import os

#Comment to insert in the backup file
comments="Change wrt f140w_FINAL_newLLM_II: same lens prior for theta_E, e1,e2, theta_E_pert and shear parameters " 
more_freedom = 3. #1. #MOD_FREEDOM -> affects only the source and perturer light parameters 

R_sersic_source = 0.057
sg_R_sersic_source = 0.0615


e1_source = 0.02
sg_e1_source= 0.1

e2_source = -0.353
sg_e2_source = 0.09

x_source = 0.466 
sg_x_source =0.016

y_source = -1.717
sg_y_source = 0.02

R_sersic_lens = 0.312
sg_R_sersic_lens = 0.04

n_sersic_lens = 3.63 
sg_n_sersic_lens= 0.4

#Modify here to indicate the image and the psf files and the model of lens light
data_path="./data/F140W/"
image_name="F140W_new_lens.fits"
err_name="e.F140W_new_lens.fits"

psf_name="small_sub_sym_psf_0001_F140W_qso_sub.fits"
err_psf="e.small_sub_sym_psf_0001_F140W_qso_sub.fits"

pssf = 1 #point source supersampling factor -> how much the PSF was supersampled
ssks = 15 #it should be the size of the non-supersampled kernel in pixel  #>note that this changed due to new psf
lens_light_model_name ="F140W_llm.fits"
#bool_bin = True #to bin the psf in the code
bool_bin=False
fixed_mag =[False] #Fixed magnification (one element per image)
sub = True #subtracted light
filter_name = "f140w"
#"cosmetics": coutout for the printed images
v_min,v_max=-1,1
e_v_min,e_v_max =-1,0
res_min, res_max= -4,4

#Coordinates
############
#CD matrix
# from degree to arcsec
CD1_1   = 1.04141997294039E-05*3600                                                  
CD1_2   = 3.40688292375962E-05*3600                                                  
CD2_1   = 3.40688292375962E-05*3600                                                  
CD2_2   = -1.0414199729403E-05*3600

from numpy import array as arr
transform_pix2angle = arr([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
ra_at_xy_0,dec_at_xy_0 = 0,0

# Parameter for the lens
center_x_pix,center_y_pix = 39.129, 44.234
dx,dy=50,50 # this must be the center pixel of the image (eg: 3 pixels [][x][] -> 2, then 1 bc python counting )
#For the perturber
x_pert_pix,y_pert_pix = 37, 29.4463
# Parameters for the source
# order A (top), C (right), B (bottom),D (left)
x_im_pix = [53.329, 35.7093, 25.3736, 41.718]
y_im_pix = [33.4973, 32.8209, 42.252, 53.859]
z_source = 2.737
z_lens= 0.407 

#Settings for cut image
edge = None

#MOD_SHIFT 
# in ra,dec considering the rotation and deformation of the reference system
mod_fr = "We consider the absolute frame of reference"
x_image, y_image, x_pert, y_pert, center_x, center_y, pix_scale, ra_at_xy_0,dec_at_xy_0 = \
	conv_5 (x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,transform_pix2angle,ra_at_xy_0,dec_at_xy_0)

#Mask coord
###########
x_mask, y_mask =[65.94,66.62,76.22],[36.41,33.67,35.8] # coord in pixel in the HST image 
#radius of the mask for the object(s) in pixels
r_mask = [2,3,2]
x_mask_conv,y_mask_conv = array(x_mask)-1,array(y_mask)-1

new_mask="correct the mask of the center"
#We implement the mask of the center
x_mask_cen,y_mask_cen =	array([40,41]),\
			array([45,44])
rad=			[3,1]


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

kwargs_lower_lens.append({'center_x': center_x-8.*pix_scale,
                          'center_y':center_y-8.*pix_scale,
                          'theta_E': min_theta_E, 
                          'e1': min_e1, 
                          'e2': min_e2})

kwargs_upper_lens.append({'center_x': center_x+8.*pix_scale, 
                          'center_y':center_y+8.*pix_scale,
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
			              'center_x': x_pert-2*pix_scale,
                          'center_y': y_pert-2*pix_scale})

kwargs_upper_lens.append({'theta_E': 0.35, 
		               	  'center_x': x_pert+2*pix_scale,
                          'center_y': y_pert+2*pix_scale})

# from f814w and f160
#Shear
fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
kwargs_lens_init.append({'gamma_ext' : 0.112, 'psi_ext': -0.41*pi})
kwargs_lens_sigma.append({'gamma_ext': 0.05, 'psi_ext': pi/50.})
kwargs_lower_lens.append({'gamma_ext': 0.00, 'psi_ext':-0.5*pi})
kwargs_upper_lens.append({'gamma_ext': 0.40, 'psi_ext':-.2*pi})
			
#Source light parameters
fixed_source = [{}]
kwargs_source_init = [{'amp' : 0,'R_sersic': R_sersic_source ,
				 'n_sersic': 4,
				 'e1': e1_source, 
				 'e2': e2_source, 
				 'center_x': x_source, 
				 'center_y': y_source}]
kwargs_source_sigma = [{'amp': 0,'R_sersic': sg_R_sersic_source*3, 
				'n_sersic': 1., 
				'e1': sg_e1_source*3 , 
				'e2': sg_e2_source*3,
				'center_x': sg_x_source*3, 
				'center_y': sg_y_source*3}]
kwargs_lower_source = [{'amp': 0,'R_sersic': R_sersic_source - (sg_R_sersic_source*3*3*more_freedom) , 
				'n_sersic': 0.5, 
				'e1': e1_source -(sg_e1_source*3*3*more_freedom), 
				'e2': e2_source - (sg_e2_source*3*3*more_freedom), 
				'center_x': x_source - (sg_x_source*3*3*more_freedom), 
				'center_y': y_source - (sg_y_source*3*3*more_freedom)}]

kwargs_upper_source = [{'amp': 0,'R_sersic': R_sersic_source + (sg_R_sersic_source*3*3*more_freedom) , 
				'n_sersic': 6., 
				'e1': e1_source + (sg_e1_source*3*3*more_freedom), 
				'e2': e2_source + (sg_e2_source*3*3*more_freedom), 
				'center_x': x_source + (sg_x_source*3*3*more_freedom), 
				'center_y': y_source + (sg_y_source*3*3*more_freedom)}]

#Ps light parameters 
fixed_ps = [{}]
lum_ps = array([0,0,0,0])
kwargs_ps_init = [{'ra_image': x_image, 
		'dec_image': y_image,
		"point_amp": lum_ps}]
bound_pix_ps=.5
kwargs_ps_sigma = [{'ra_image' : 1*pix_scale*ones_like(x_image), 
                    'dec_image': 1*pix_scale*ones_like(x_image),
                    "point_amp":lum_ps  }]
kwargs_lower_ps = [{'ra_image' : x_image-bound_pix_ps*pix_scale,
                    'dec_image': y_image-bound_pix_ps*pix_scale,
                    "point_amp": lum_ps}]
kwargs_upper_ps = [{'ra_image' : x_image+bound_pix_ps*pix_scale,
                    'dec_image': y_image+bound_pix_ps*pix_scale,
                    "point_amp": lum_ps}]
#Lens light parameters
fixed_lens_light = []
kwargs_lens_light_init = []
kwargs_lens_light_sigma = []
kwargs_lower_lens_light = []
kwargs_upper_lens_light = []
# Perturber light profile
bound_pix_pert=5
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					'n_sersic':n_sersic_lens , 
					'center_x':x_pert,
					'center_y':y_pert})
kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					'n_sersic': sg_n_sersic_lens*3, 
					'center_x':6.*pix_scale,
					'center_y':6.*pix_scale}) 
kwargs_lower_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens - (sg_R_sersic_lens*3*3*more_freedom), 
					'n_sersic': 0.,
	                                'center_x':x_pert - bound_pix_pert*pix_scale,
	                                'center_y':y_pert - bound_pix_pert*pix_scale}) 
kwargs_upper_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens + (sg_R_sersic_lens*3*3*more_freedom), 
					'n_sersic': n_sersic_lens + (sg_n_sersic_lens*3*more_freedom),
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
	from produce_kwargs_result import standard
	return standard(samples_mcmc,param_mcmc,i)
	

def str_param(param):
    from get_param_title import standard
    return standard(param)
