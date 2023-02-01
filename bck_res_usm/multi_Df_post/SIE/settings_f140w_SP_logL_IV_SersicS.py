#Initial data parameter settings for f140w with lens light subtracted
#####################################################################
import warnings
import numpy as np
from conversion import conv_5 as conv
from conversion import find_barycenter
from same_prior_ellipt_logL_II import lens_prior
from image_manipulation import get_transf_matrix

class setting():

    warnings.simplefilter("ignore")

    def __init__(self):
        #Comment to insert in the backup file
        self.comments="Change wrt settings_f140w_SP_ellipt_logL_V.py: use sersic source"
                
        # copy from f160w_logL_III
        R_sersic_source = 0.298
        sg_R_sersic_source = 0.02

        n_sersic_source = 1.35
        sg_n_sersic_source= .2
        
        x_source = 0.48
        y_source = -1.71
        sg_x_source = 0.02
        sg_y_source = 0.02       
        
        
        R_sersic_lens    = 0.312
        sg_R_sersic_lens = 0.04
        
        n_sersic_lens    = 3.63 
        sg_n_sersic_lens = 0.4

        #Modify here to indicate the image and the psf files and the model of lens light
        self.data_path  = "./data/F140W/"
        self.image_name = "F140W_new_lens.fits"
        self.err_name   = "e.F140W_new_lens.fits"
        
        self.psf_name   = "small_sub_sym_psf_0001_F140W_qso_sub.fits"
        self.err_psf    = "e.small_sub_sym_psf_0001_F140W_qso_sub.fits"
        
        self.pssf = 1 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 15 #it should be the size of the non-supersampled kernel in pixel  #>note that this changed due to new psf
        self.lens_light_model_name ="F140W_llm.fits"
        self.fixed_mag = [False] #Fixed magnification (one element per image)
        self.sub = False #subtracted light
        self.already_sub = False
        self.filter_name = "f140w"
        #"cosmetics": coutout for the printed images
        self.v_min  , self.v_max   = -1,1
        self.e_v_min, self.e_v_max = -1,0
        self.res_min, self.res_max = -4,4

        #Coordinates
        ############
        #CD matrix
        self.transform_pix2angle           = get_transf_matrix(self.data_path+self.image_name,in_arcsec=True)
        self.ra_at_xy_0,self.dec_at_xy_0   = 0,0

        # Parameter for the lens
        center_x_pix,center_y_pix = 40.1009, 45.2448
        #For the perturber
        x_pert_pix,y_pert_pix     = 37, 29.4463
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [53.329, 35.7093, 25.3736, 41.718]
        y_im_pix = [33.4973, 32.8209, 42.252, 53.859]
        self.z_source = 2.737
        self.z_lens   = 0.407 

        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0, self.dec_at_xy_0 = \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,self.ra_at_xy_0,self.dec_at_xy_0)

        #About this model
        self.WS = False #without source modelled
        self.ellipt=False#no elliptical source
        self.CP = False #change Profile

        #Mask coord
        ###########
        x_mask, y_mask          = [65.94,66.62,76.22],[36.41,33.67,35.8] # coord in pixel in the HST image 
        #radius of the mask for the object(s) in pixels
        self.r_mask             =  [2,3,2]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        #We implement the mask of the center
        x_mask_cen,y_mask_cen           = np.array([40,41]),np.array([45,44])
        self.x_mask_cen,self.y_mask_cen = x_mask_cen-1,y_mask_cen-1
        self.rad_mask_cen               = [3,1]

        self.R_tot_mask                 = 3./self.pix_scale # 3 "-> pixel
        self.new_center_mask_out        = find_barycenter(x_im_pix,y_im_pix)

        #Model parameters
        #################
        #Lens mass profile
        bound_pix_mainlens = 2.
        bound_mainlens     = bound_pix_mainlens*self.pix_scale

        fixed_lens        = []
        kwargs_lens_init  = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []

        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E': lens_prior.theta_E,
                                 'center_x': self.center_x, 
			                     'center_y': self.center_y,
			                     'e1': lens_prior.e1,
		         	             'e2': lens_prior.e2})
        kwargs_lens_sigma.append({'center_x': .9*self.pix_scale, 
			                      'center_y': .9*self.pix_scale,
                                  'theta_E': lens_prior.sig_theta_E, 
                                  'e1': lens_prior.sig_e1, 
                                  'e2': lens_prior.sig_e2})

        kwargs_lower_lens.append({'center_x': self.center_x-bound_mainlens,
                                  'center_y': self.center_y-bound_mainlens,
                                  'theta_E': lens_prior.min_theta_E, 
                                  'e1': lens_prior.min_e1, 
                                  'e2': lens_prior.min_e2})

        kwargs_upper_lens.append({'center_x': self.center_x+bound_mainlens, 
                                  'center_y': self.center_y+bound_mainlens,
                                  'theta_E': lens_prior.max_theta_E, 
                                  'e1': lens_prior.max_e1, 
                                  'e2': lens_prior.max_e2})
                                                             
        # Perturber mass profile
        bound_pix_pert = 2.
        bound_pert     = bound_pix_pert*self.pix_scale
        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E' : lens_prior.theta_pert ,
                                 'center_x': self.x_pert, 
                                 'center_y': self.y_pert})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_pert , 
			                      'center_x': 6.*self.pix_scale,
                                  'center_y': 6.*self.pix_scale})

        kwargs_lower_lens.append({'theta_E': lens_prior.min_theta_pert,
			                      'center_x': self.x_pert-bound_pert,
                                  'center_y': self.y_pert-bound_pert})

        kwargs_upper_lens.append({'theta_E': lens_prior.max_theta_pert, 
		                       	  'center_x': self.x_pert+bound_pert,
                                  'center_y': self.y_pert+bound_pert})

        #Shear
        fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
        kwargs_lens_init.append({'gamma_ext' : lens_prior.gamma_ext,     'psi_ext': lens_prior.psi_ext})
        kwargs_lens_sigma.append({'gamma_ext': lens_prior.sig_gamma_ext, 'psi_ext': lens_prior.sig_psi_ext})
        kwargs_lower_lens.append({'gamma_ext': lens_prior.min_gamma_ext, 'psi_ext': lens_prior.min_psi_ext})
        kwargs_upper_lens.append({'gamma_ext': lens_prior.max_gamma_ext, 'psi_ext': lens_prior.max_psi_ext})
                          
	        
        #Source light parameters
        # copy from f160w_SP_logL_III
        fixed_source = [{}]
        kwargs_source_init = [{'amp' : 0,'R_sersic': R_sersic_source ,
				                'n_sersic': n_sersic_source, 
				                'center_x': x_source, 
				                'center_y': y_source}]
        kwargs_source_sigma = [{'amp': 0,'R_sersic': sg_R_sersic_source*3, 
				                'n_sersic': n_sersic_source*3, 
				                'center_x': sg_x_source*3, 
				                'center_y': sg_y_source*3}]
        kwargs_lower_source = [{'amp': 0,'R_sersic': R_sersic_source - (sg_R_sersic_source*9) , 
				                'n_sersic': 1., 
				                'center_x': x_source - (sg_x_source*9), 
				                'center_y': y_source - (sg_y_source*9)}]
        kwargs_upper_source = [{'amp': 0,'R_sersic': R_sersic_source + (sg_R_sersic_source*9) , 
				                'n_sersic': 7., 
				                'center_x': x_source + (sg_x_source*9), 
				                'center_y': y_source + (sg_y_source*9)}]


        #Ps light parameters 
        fixed_ps = [{}]
        lum_ps = np.array([0,0,0,0])
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image': self.y_image,"point_amp": lum_ps}]
        bound_pix_ps = 1.5
        bound_ps     = bound_pix_ps*self.pix_scale
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),
                            "point_amp":lum_ps}]
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_ps,
                            'dec_image': self.y_image-bound_ps,
                            "point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : self.x_image+bound_ps,
                            'dec_image': self.y_image+bound_ps,
                            "point_amp": lum_ps}]
        
        # Lens light parameters
        fixed_lens_light = [{}] 
        kwargs_lens_light_init =[ {'amp': 0.,'R_sersic':1.46,
				        'n_sersic':4.5, 
		         		'e1': 0.26 ,
		         		'e2': -0.07,
		         		'center_x': self.center_x, 
		         		'center_y': self.center_y}]
        kwargs_lens_light_sigma=[{'amp': 0.,'R_sersic':.2, 
				        'n_sersic':0.16, 
				        'e1': 0.026*3,
				        'e2': 0.015*3,
				        'center_x': .8*self.pix_scale, 
				        'center_y': .8*self.pix_scale}]
        kwargs_lower_lens_light=[{'amp': 0.,'R_sersic': 1,
				        'n_sersic':2.7,
			         	'e1': 0.213-(0.026*9) ,
			         	'e2': -0.089-( 0.015*9),
				        'center_x': self.center_x-bound_mainlens,
				        'center_y': self.center_y-bound_mainlens}]
        kwargs_upper_lens_light = [{'amp':0.,'R_sersic':2, 
				        'n_sersic':6,
			         	'e1': 0.213+(0.026*9) ,
			         	'e2': -0.089+( 0.015*9),
				        'center_x': self.center_x+bound_mainlens,
				        'center_y': self.center_y+bound_mainlens}]
				        
        # Perturber light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
    					               'n_sersic':n_sersic_lens, 
					                   'center_x':self.x_pert,
					                   'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					                    'n_sersic': sg_n_sersic_lens*3, 
					                    'center_x':6.*self.pix_scale,
					                    'center_y':6.*self.pix_scale}) 
        kwargs_lower_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens - (sg_R_sersic_lens*9), 
					                    'n_sersic': 0.,
	                                    'center_x':self.x_pert - bound_pert,
	                                    'center_y':self.y_pert - bound_pert}) 
        kwargs_upper_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens + (sg_R_sersic_lens*9), 
            					        'n_sersic':n_sersic_lens + (sg_n_sersic_lens*3),
                                        'center_x':self.x_pert + bound_pert,
	                                    'center_y':self.y_pert + bound_pert}) 

        #Uniform light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp':1.5})
        kwargs_lens_light_sigma.append({'amp':3})
        kwargs_lower_lens_light.append({'amp':0.})
        kwargs_upper_lens_light.append({'amp':10.}) 

        lens_params       = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        ps_params         = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]
        source_params     = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]


        # Model parameters
        self.lens_params       = lens_params
        self.ps_params         = ps_params
        self.lens_light_params = lens_light_params
        self.source_params     = source_params


    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
        from produce_kwargs_result import with_main_lens_light as produce_kwargs_result
        return produce_kwargs_result(samples_mcmc,param_mcmc,i)

    def str_param(self,param):
        from get_param_title import with_main_lens_light as str_param
        return str_param(param)
