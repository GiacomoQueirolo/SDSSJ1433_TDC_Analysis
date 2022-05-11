#Initial data parameter settings for f814w with lens light subtracted
#####################################################################
from numpy import array, inf, nan, pi,ones_like
from conversion import conv, conv_5
import os

#Comment to insert in the backup file
comments="Change wrt f814w_final_VI_ws:  same lens prior for theta_E, e1,e2, theta_E_pert and shear parameters " 
#Modify here to indicate the image and the psf files and the model of lens light
data_path="./data/F814W/"
image_name="F814W_qso.fits"
err_name="e.F814W_qso.fits"
#psf_name="psf_0005_hst_f814w.fits"
psf_name="psf_0005_F814W.fits"
err_psf="e.psf_0001_F814W.fits"
pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
ssks = 41#it should be the size of the non-supersampled kernel in pixel
lens_light_model_name=None
bool_bin=False
fixed_mag =[False] #Fixed magnification (one element per image)
sub = True #subtracted light
filter_name = "f814w"
#"cosmetics": coutout for the printed images
v_min,v_max=-4,1
e_v_min,e_v_max =-3,0
res_min, res_max= -3,3

#Coordinates
############
#CD matrix
# from degree to arcsec

CD1_1   = -8.1353962040008E-06*3600                                               
CD1_2   = 7.41198940931543E-06*3600                                               
CD2_1   = 7.41198940931543E-06*3600                                             
CD2_2   = 8.13539620400086E-06*3600
 
from numpy import array as arr
transform_pix2angle = arr([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
ra_at_xy_0,dec_at_xy_0 = 0,0

# Parameter for the lens
center_x_pix,center_y_pix = 62.6044,87.1731 

theta_E = 1.645 
#For the perturber
x_pert_pix,y_pert_pix = 104.539,56
# Parameters for the source
# order A (top), C (right), B (bottom),D (left)
x_im_pix = [114.882,39.794,51.0978,92.899]
y_im_pix = [109.799,103.96,39.6647,57.13]
z_source = 2.737
z_lens= 0.407
 
#Settings for cut image
edge = None

#Conversion in angle in the cut image

#MOD_SHIFT 
# in ra,dec considering the rotation and deformation of the reference system
mod_fr = "We consider the absolute frame of reference"
x_image, y_image, x_pert, y_pert, center_x, center_y, pix_scale, ra_at_xy_0,dec_at_xy_0= \
	conv_5 (x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,transform_pix2angle,ra_at_xy_0,dec_at_xy_0)
#Mask coord
###########
x_mask, y_mask =[118],[148] # coord in pixel in the HST image 
#radius of the mask for the object in pixels
r_mask = [20.824865]
x_mask_conv,y_mask_conv = array(x_mask)-1,array(y_mask)-1

new_mask="we consider a new masking method"
#We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
x_mask_cen,y_mask_cen =	array([0,-1,+1,0,+1,0,-1])+center_x_pix,\
			array([0,+1,-1,-2,-2,-1,0])+center_y_pix
rad=			[2,1,1,1,1,1,3]

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

kwargs_lens_sigma.append({'center_x': .8*pix_scale, 
			            'center_y':.8*pix_scale,
                        'theta_E': sig_theta_E, 
                        'e1': sig_e1, 
                        'e2': sig_e2})

kwargs_lower_lens.append({'center_x': center_x-4.*pix_scale,
                          'center_y':center_y-4.*pix_scale,
                          'theta_E': min_theta_E, 
                          'e1': min_e1, 
                          'e2': min_e2})

kwargs_upper_lens.append({'center_x': center_x+6.*pix_scale, 
                          'center_y':center_y+3.*pix_scale,
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
			              'center_x': 2.*pix_scale,
                          'center_y': 2.*pix_scale})

kwargs_lower_lens.append({'theta_E': 0.,
			              'center_x': x_pert-6*pix_scale,
                          'center_y': y_pert-6*pix_scale})

kwargs_upper_lens.append({'theta_E': 0.35, 
		               	  'center_x': x_pert+6*pix_scale,
                          'center_y': y_pert+6*pix_scale})

# from f814w and f160
#Shear
fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
kwargs_lens_init.append({'gamma_ext' : 0.112, 'psi_ext': -0.41*pi})
kwargs_lens_sigma.append({'gamma_ext': 0.05, 'psi_ext': pi/50.})
kwargs_lower_lens.append({'gamma_ext': 0.00, 'psi_ext':-0.5*pi})
kwargs_upper_lens.append({'gamma_ext': 0.40, 'psi_ext':-.2*pi})


#Ps light parameters
fixed_ps = [{}]
#to consider the amplitude of the point source, we measure the maximum count in the image for each
lum_ps = array([24.,22.,30.,5.])
#Lower and upper boundary position freedom in pixel wrt the input value
bound_pix_ps=.6 #1
kwargs_ps_init = [{'ra_image': x_image, 'dec_image': y_image,"point_amp": lum_ps}]
kwargs_ps_sigma = [{'ra_image': .2*pix_scale*ones_like(x_image), 
                    'dec_image': .2*pix_scale*ones_like(x_image),"point_amp":lum_ps/10. }]
kwargs_lower_ps = [{'ra_image': x_image-bound_pix_ps*pix_scale, 
                    'dec_image': y_image-bound_pix_ps*pix_scale,"point_amp": .1*lum_ps}]
kwargs_upper_ps = [{'ra_image': x_image+bound_pix_ps*pix_scale, 
                    'dec_image': y_image+bound_pix_ps*pix_scale,"point_amp": 10*lum_ps}]

#Lens light parameters
fixed_lens_light = []
kwargs_lens_light_init = []
kwargs_lens_light_sigma = []
kwargs_lower_lens_light = []
kwargs_upper_lens_light = []
# Perturber light profile
#Lower and upper boundary position freedom in pixel wrt the input value
bound_pix_pert=7 #5
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp':14.,'R_sersic':0.375, 
				'n_sersic':4.09, 
                               'center_x':x_pert,'center_y':y_pert})
kwargs_lens_light_sigma.append({'amp':10.,'R_sersic':.1,  
				'n_sersic':.6, 
                                'center_x':4.*pix_scale,'center_y':4.*pix_scale})
kwargs_lower_lens_light.append({'amp':3,'R_sersic':0.1, 
				'n_sersic':1.,
                                'center_x':x_pert-bound_pix_pert*pix_scale,
                                'center_y':y_pert-bound_pix_pert*pix_scale})
kwargs_upper_lens_light.append({'amp':25.,'R_sersic':1.,
				 'n_sersic':6,
				  'center_x':x_pert+bound_pix_pert*pix_scale,
				  'center_y':y_pert+bound_pix_pert*pix_scale})

#Uniform light profile
fixed_lens_light.append({})
kwargs_lens_light_init.append({'amp':0.5})
kwargs_lens_light_sigma.append({'amp':2})
kwargs_lower_lens_light.append({'amp':0.})
kwargs_upper_lens_light.append({'amp':5.})




#MOD_PROD
def produce_kwargs_result(samples_mcmc,param_mcmc,i):
	from produce_kwargs_result import standard
	return standard(samples_mcmc,param_mcmc,i)
	

def str_param(param):
    from get_param_title import standard
    return standard(param)
