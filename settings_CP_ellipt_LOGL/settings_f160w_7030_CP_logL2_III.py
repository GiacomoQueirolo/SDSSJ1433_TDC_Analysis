#Initial data parameter settings for f160w 
##########################################
import warnings
import numpy  as np
from conversion import conv_5 as conv
from conversion import find_barycenter
from same_prior_ellipt_logL_II import lens_PEMD_prior as lens_prior
from image_manipulation import get_transf_matrix

class setting():

    warnings.simplefilter("ignore")

    def __init__(self):
        #Comment to insert in the backup file
        self.comments = "Change wrt settings_f160w_7030_CP_logL_II.py: correct PSF model"

        self.data_path  = "./data/F160W_7030/"
        self.image_name = "F160W_qso_7030.fits"
        self.err_name   = "e.F160W_qso_7030.fits"

        self.psf_name   = "psf_0005_maybe_7030.fits"
        self.err_psf    = "e.psf_0001_maybe_7030.fits"

        #Perturber
        R_sersic_lens    = 0.363
        sg_R_sersic_lens = 0.1
        n_sersic_lens    = 4.0 
        sg_n_sersic_lens = 0.5

        x_source = 0.48
        y_source = -1.71
        sg_x_source = 0.02
        sg_y_source = 0.02

        #from f160w_FINAL_mask_V:
        ######################################
        R_sersic_source = 0.298
        sg_R_sersic_source = 0.02

        n_sersic_source = 1.35
        sg_n_sersic_source= .2

        e1_source = 0.09
        sg_e1_source= 0.02

        e2_source = -0.1
        sg_e2_source = 0.02
        ########################################
        # Parameter for the lens
        center_x_pix,center_y_pix = 85.2789, 93.9701
        #For the perturber
        x_pert_pix,y_pert_pix = 98.0000, 84.4429
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [101.459,94.5911,81.6961,78.3847]
        y_im_pix = [100.94,84.6603,79.3016,98.9967]
        
        self.z_source = 2.737
        self.z_lens   = 0.407  

        self.pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 71	#it should be the size of the non-supersampled kernel in pixel
        #lens_light_model_name ="model_f814w_large.fits"
        #lens_light_model_name="F160W_new_llm_II.fits"
        self.lens_light_model_name= None
        #bool_bin = True #to bin the psf in the code
        self.bool_bin=False
        self.fixed_mag =[False] #Fixed magnification (one element per image)
        self.sub =False #subtracted light
        self.already_sub = False
        self.filter_name = "f160w_7030"
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max     = -1,1
        self.e_v_min,self.e_v_max = -3,0
        self.res_min,self.res_max = -6,6
        #Coordinates
        ############
        #CD matrix
        # from degree to arcsec
        #CD matrix
        self.transform_pix2angle         = get_transf_matrix(self.data_path+self.image_name,in_arcsec=True)
        self.ra_at_xy_0,self.dec_at_xy_0 = 0,0

        #Conversion in angle
        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0,self.dec_at_xy_0 = \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,self.ra_at_xy_0,self.dec_at_xy_0)
	        
        #About this model
        self.WS = False #without source modelled
        self.ellipt=False # w/o source ellipticity
        self.CP = True #change Profile

        #Mask coord
        ###########
        x_mask, y_mask =[],[] # coord in pixel in the HST image 
        #radius of the mask for the object in pixels
        self.r_mask = []
        self.x_mask,self.y_mask = np.array([]),np.array([])

        #We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
        x_mask_cen,y_mask_cen           = np.array([]),np.array([])
        self.x_mask_cen,self.y_mask_cen = x_mask_cen-1,y_mask_cen-1
        self.rad_mask_cen               = []
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
			                     'center_y':self.center_y,
			                     'gamma': lens_prior.gamma,
			                     'e1': lens_prior.e1,
		         	             'e2': lens_prior.e2})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_E, 
                                  'center_x': .9*self.pix_scale, 
			                      'center_y': .9*self.pix_scale,
			                      'gamma': lens_prior.sig_gamma,
                                  'e1': lens_prior.sig_e1, 
                                  'e2': lens_prior.sig_e2})

        kwargs_lower_lens.append({'theta_E': lens_prior.min_theta_E, 
                                  'center_x': self.center_x-bound_mainlens,
                                  'center_y': self.center_y-bound_mainlens,
			                      'gamma': lens_prior.min_gamma,
                                  'e1': lens_prior.min_e1, 
                                  'e2': lens_prior.min_e2})

        kwargs_upper_lens.append({'theta_E': lens_prior.max_theta_E, 
			                      'center_x': self.center_x+bound_mainlens, 
                                  'center_y': self.center_y+bound_mainlens,
                                  'gamma': lens_prior.max_gamma,
                                  'e1': lens_prior.max_e1, 
                                  'e2': lens_prior.max_e2})

        # Perturber mass profile
        bound_pix_pert = 2.
        bound_pert     = bound_pix_pert*self.pix_scale
        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E' : lens_prior.theta_pert,
                                 'center_x': self.x_pert, 
                                 'center_y': self.y_pert})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_pert, 
			                    'center_x' : 1.*self.pix_scale,
                                'center_y' : 1.*self.pix_scale})

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
        lum_ps = np.array([0,0,0,0])
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image': self.y_image,"point_amp": lum_ps}]
        bound_pix_ps = 1.3
        bound_ps     = bound_pix_ps*self.pix_scale
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),"point_amp":lum_ps}]
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_ps, 
                            'dec_image': self.y_image-bound_ps,"point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image': self.x_image+bound_ps, 
                            'dec_image': self.y_image+bound_ps,"point_amp": lum_ps}]
        #Lens light parameters
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
				        'center_y':self.center_y-bound_mainlens}]
        kwargs_upper_lens_light = [{'amp':0.,'R_sersic':2, 
				        'n_sersic':6,
			         	'e1': 0.213+(0.026*9) ,
			         	'e2': -0.089+( 0.015*9),
				        'center_x': self.center_x+bound_mainlens,
				        'center_y': self.center_y+bound_mainlens}]


        # Perturber light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					                    'n_sersic':n_sersic_lens , 
					                    'center_x':self.x_pert,
				                    	'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					                    'n_sersic': sg_n_sersic_lens*3, 
					                    'center_x':1.*self.pix_scale,
					                    'center_y':1.*self.pix_scale}) 
        kwargs_lower_lens_light.append({'amp':0.,'R_sersic':.1, 
			                        	'n_sersic': 1,
                                        'center_x':self.x_pert - bound_pert,
                                        'center_y':self.y_pert - bound_pert}) 
        kwargs_upper_lens_light.append({'amp':0.,'R_sersic':2, 
				                        'n_sersic': 5,
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
        
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
	    from produce_kwargs_result import newProfile_with_main_lens_light as produce_kwargs_result
	    return produce_kwargs_result(samples_mcmc,param_mcmc,i)
	    
    def str_param(self,param):
        from get_param_title import newProfile_with_main_lens_light as str_param
        return str_param(param)
