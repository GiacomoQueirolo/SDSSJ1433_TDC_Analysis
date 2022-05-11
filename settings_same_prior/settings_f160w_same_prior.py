#Initial data parameter settings for f160w 
##########################################
from numpy import array, inf, nan, pi,ones_like
from conversion import conv, conv_5
import os

#Comment to insert in the backup file
comments="Change wrt f160w_FINAL_mask_V:same lens prior for theta_E, e1,e2, theta_E_pert and shear parameters "   

R_sersic_source = 0.057
sg_R_sersic_source = 0.0615

n_sersic_source = 4
sg_n_sersic_source= 1

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
data_path="./data/F160W/"
image_name="qso_F160W_arno.fits"
err_name="e.qso_F160W_arno.fits"

psf_name="psf_0005_F160W_3stars.fits"
err_psf="e.psf_0001_F160W_3stars.fits"

pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
ssks = 31	#it should be the size of the non-supersampled kernel in pixel
#lens_light_model_name ="model_f814w_large.fits"
lens_light_model_name= None
#bool_bin = True #to bin the psf in the code
bool_bin=False
fixed_mag =[False] #Fixed magnification (one element per image)
sub =False #subtracted light
filter_name = "f160w"
#"cosmetics": coutout for the printed images
v_min,v_max=-1,1
e_v_min,e_v_max =-3,0
res_min, res_max= -6,6

#Coordinates
############
#CD matrix
# from degree to arcsec

CD1_1   = -2.6334282913251E-05*3600                                                 
CD1_2   = 2.39926285568969E-05*3600                                           
CD2_1   = 2.39926285568969E-05*3600                                            
CD2_2   = 2.63342829132518E-05*3600 
 
from numpy import array as arr
transform_pix2angle = arr([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
ra_at_xy_0,dec_at_xy_0 = 0,0

# Parameter for the lens
center_x_pix,center_y_pix = 52.1178,52.2733
#For the perturber
x_pert_pix,y_pert_pix = 65.2,42.6
# Parameters for the source
# order A (top), C (right), B (bottom),D (left)
x_im_pix = [68.2652, 62.2084, 48.8707, 45.1079]
y_im_pix = [59.1228, 42.843, 37.5594, 57.4156]

z_source = 2.737
z_lens= 0.407

#Settings for cut image
edge = None	

#Conversion in angle in the cut image

#MOD_SHIFT 
# in ra,dec considering the rotation and deformation of the reference system
mod_fr = "We consider the absolute frame of reference"
x_image, y_image, x_pert, y_pert, center_x, center_y, pix_scale, ra_at_xy_0,dec_at_xy_0 = \
	conv_5(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,transform_pix2angle,ra_at_xy_0,dec_at_xy_0)

#Mask coord
###########
x_mask, y_mask =[72.5,77.5,82.8,76.1,88.4,47,54,69.,78,101,52],[71.8,69.87,68.2,80.48,46.4,100,101,32,26,52,52] # coord in pixel in the HST image 
#radius of the mask for the object in pixels
r_mask = [3.9,2.33,2.33,3.42,3.42,3,2,2.5,2,1,3]
x_mask_conv,y_mask_conv = array(x_mask)-1,array(y_mask)-1

new_mask="we consider a new masking method"
#We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
x_mask_cen,y_mask_cen =	array([]),\
			array([])
rad=			[]

#Model parameters

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

kwargs_lens_sigma.append({'center_x': .6*pix_scale, 
			            'center_y':.6*pix_scale,
                        'theta_E': sig_theta_E, 
                        'e1': sig_e1, 
                        'e2': sig_e2})

kwargs_lower_lens.append({'center_x': center_x-2.*pix_scale,
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
			              'center_x': 1.*pix_scale,
                          'center_y': 1.*pix_scale})

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
kwargs_source_init = [{'amp' : 0,'R_sersic': R_sersic_source ,
				 'n_sersic': n_sersic_source, 
				 'e1': e1_source, 
				 'e2': e2_source, 
				 'center_x': x_source, 
				 'center_y': y_source}]
kwargs_source_sigma = [{'amp': 0,'R_sersic': sg_R_sersic_source*3, 
				'n_sersic': n_sersic_source*3, 
				'e1': sg_e1_source*3 , 
				'e2': sg_e2_source*3,
				'center_x': sg_x_source*3, 
				'center_y': sg_y_source*3}]
kwargs_lower_source = [{'amp': 0,'R_sersic': R_sersic_source - (sg_R_sersic_source*9) , 
				'n_sersic': 1., 
				'e1': e1_source -.9, 
				'e2': e2_source - .9, 
				'center_x': x_source - (sg_x_source*9), 
				'center_y': y_source - (sg_y_source*9)}]

kwargs_upper_source = [{'amp': 0,'R_sersic': R_sersic_source + (sg_R_sersic_source*9) , 
				'n_sersic': 7., 
				'e1': e1_source + (sg_e1_source*9), 
				'e2': e2_source + (sg_e2_source*9), 
				'center_x': x_source + (sg_x_source*9), 
				'center_y': y_source + (sg_y_source*9)}]

#Ps light parameters 
fixed_ps = [{}]
lum_ps = array([0,0,0,0])
kwargs_ps_init = [{'ra_image': x_image, 
		'dec_image': y_image,
		"point_amp": lum_ps}]
bound_pix_ps=1.3
kwargs_ps_sigma = [{'ra_image' : 1*pix_scale*ones_like(x_image), 
                    'dec_image': 1*pix_scale*ones_like(x_image),"point_amp":lum_ps  }]
kwargs_lower_ps = [{'ra_image' : x_image-bound_pix_ps*pix_scale, 
                    'dec_image': y_image-bound_pix_ps*pix_scale,"point_amp": lum_ps}]
kwargs_upper_ps = [{'ra_image' : x_image+bound_pix_ps*pix_scale, 
                    'dec_image': y_image+bound_pix_ps*pix_scale,"point_amp": lum_ps}]
#Lens light parameters
fixed_lens_light = [{}] 
kwargs_lens_light_init =[ {'amp': 0.,'R_sersic':1.46,
				'n_sersic':4.5, 
		 		'e1': 0.26 ,
		 		'e2': -0.07,
		 		'center_x': center_x, 
		 		'center_y': center_y}]
kwargs_lens_light_sigma=[{'amp': 0.,'R_sersic':.2, 
				'n_sersic':0.16, 
				'e1': 0.026*3,
				'e2': 0.015*3,
				'center_x': .8*pix_scale, 
				'center_y': .8*pix_scale}]
kwargs_lower_lens_light=[{'amp': 0.,'R_sersic': 1,
				'n_sersic':2.7,
			 	'e1': 0.213-(0.026*9) ,
			 	'e2': -0.089-( 0.015*9),
				'center_x': center_x-2.*pix_scale,
				'center_y':center_y-2.*pix_scale,}]
kwargs_upper_lens_light = [{'amp':0.,'R_sersic':2, 
				'n_sersic':6,
			 	'e1': 0.213+(0.026*9) ,
			 	'e2': -0.089+( 0.015*9),
				'center_x': center_x+2.*pix_scale,
				'center_y':center_y+2.*pix_scale}]
# Perturber light profile
bound_pix_pert=3
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					'n_sersic':n_sersic_lens , 
					'center_x':x_pert,
					'center_y':y_pert})
kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					'n_sersic': sg_n_sersic_lens*3, 
					'center_x':1.*pix_scale,
					'center_y':1.*pix_scale}) 
kwargs_lower_lens_light.append({'amp':0.,'R_sersic':.1, 
					'n_sersic': 1,
	                                'center_x':x_pert - bound_pix_pert*pix_scale,
	                                'center_y':y_pert - bound_pix_pert*pix_scale}) 
kwargs_upper_lens_light.append({'amp':0.,'R_sersic':2, 
					'n_sersic': 5,
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
