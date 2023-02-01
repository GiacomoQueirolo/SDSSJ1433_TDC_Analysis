#Initial data parameter settings for f475w with lens light subtracted
#####################################################################
import warnings
import numpy as np # array, pi,ones_like
from conversion import conv_5 as conv
from conversion import find_barycenter
from same_prior_ellipt_logL_II import lens_PEMD_prior as lens_prior
from image_manipulation import extract_phi_ll,extract_q_ll,get_transf_matrix

class setting(): 

    warnings.simplefilter("ignore")
   
    def __init__(self):
        
        #Comment to insert in the backup file
        self.comments = "Change wrt settings_f475w_CP_logL_II_ws: mod for MOD_CUSTOM_LIKE_II"

        R_sersic_lens    = 0.363
        sg_R_sersic_lens = 0.1

        n_sersic_lens    = 4.0 
        sg_n_sersic_lens = 0.5

        x_source    = 0.48
        y_source    = -1.71
        sg_x_source = 0.02
        sg_y_source = 0.02


        #Modify here to indicate the image and the psf files and the model of lens light
        
        self.data_path   = "./data/F475X/llm_f475w_220527/"
        self.image_name  = "f475w_qso_220527.fits"
        self.err_name    = "e."+self.image_name
        self.psf_name    = "psf_0005_F475X_2qso.fits"
        self.err_psf     = "e.psf_0001_F475X_2qso.fits"
        self.filter_name = "f475w"
        self.pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 41 #it should be the size of the non-supersampled kernel in pixel
        self.bool_bin  = False
        self.fixed_mag = [False] #Fixed magnification (one element per image)
        
        self.sub = True #subtracted light
        self.lens_light_model_name =  'f475w_llm_220527.fits' 
        self.lens_light_model_name_full = "model_tot_f475w_220527.fits" 
        self.already_sub = True # is the lens light model already subtracted?

        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max      = -4,1
        self.e_v_min,self.e_v_max  = -3,0
        self.res_min, self.res_max = -6,6
        #Coordinates
        ############
        #CD matrix
        self.transform_pix2angle         = get_transf_matrix(self.data_path+self.image_name,in_arcsec=True)
        self.ra_at_xy_0,self.dec_at_xy_0 = 0,0

        # Parameter for the lens
        center_x_pix,center_y_pix = 71.6775,91.8765

        #For the perturber
        x_pert_pix,y_pert_pix = 113.,60.5125
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [123.958, 101.961, 60.2043, 48.8816]
        y_im_pix = [114.5,61.828,44.3867,108.661]
        
        self.z_source = 2.737
        self.z_lens   = 0.407
        #Conversion in angle
        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0,self.dec_at_xy_0 = \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,self.ra_at_xy_0,self.dec_at_xy_0)

        #About this model
        self.WS = True #without source modelled
        self.CP = True #change Profile

        # MOD_CUSTOM_LIKE_II
        self.phi_ll = extract_phi_ll(model_name=self.lens_light_model_name_full,setting=self,min_ra=lens_prior.theta_E)
        self.q_ll   = extract_q_ll(model_name=self.lens_light_model_name_full,setting=self,min_ra=lens_prior.theta_E)

        #Mask coord
        ###########
        x_mask, y_mask = [137],[155] # coord in pixel in the HST image 
        #radius of the mask for the object in pixels
        self.r_mask = [10]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        #We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
        self.x_mask_cen,self.y_mask_cen = np.array([72])-1,np.array([92])-1
        self.rad_mask_cen               = [4]
        self.R_tot_mask                 = 3./self.pix_scale # 3 "-> pixel
        self.new_center_mask_out        = find_barycenter(x_im_pix,y_im_pix)

        #Model parameters
        #################

        #Lens mass profile
        bound_pix_mainlens = 4.
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
        bound_pix_pert = 4.
        bound_pert     = bound_pix_pert*self.pix_scale

        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E' : lens_prior.theta_pert,
                                 'center_x': self.x_pert, 
                                 'center_y': self.y_pert})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_pert, 
			                      'center_x': 2.*self.pix_scale,
                                  'center_y': 2.*self.pix_scale})

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
        # NONE
        #Ps light parameters 
        fixed_ps = [{}]
        lum_ps   = np.array([0,0,0,0])
        bound_pix_ps = 1.8 #.6
        bound_ps     = bound_pix_ps*self.pix_scale
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image':self.y_image,"point_amp": lum_ps}]
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),"point_amp": lum_ps}]
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_ps, 'dec_image': self.y_image-bound_ps,"point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : self.x_image+bound_ps, 'dec_image': self.y_image+bound_ps,"point_amp": lum_ps}]
        #Lens light parameters
        fixed_lens_light        = []
        kwargs_lens_light_init  = []
        kwargs_lens_light_sigma = []
        kwargs_lower_lens_light = []
        kwargs_upper_lens_light = []
        # Perturber light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					                    'n_sersic':n_sersic_lens , 
					                    'center_x':self.x_pert,
					                    'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					                    'n_sersic': sg_n_sersic_lens*3, 
					                    'center_x':6.*self.pix_scale,
					                    'center_y':6.*self.pix_scale}) 
        kwargs_lower_lens_light.append({'amp':0.,'R_sersic': 0, 
			                        	'n_sersic': 0,
                                        'center_x':self.x_pert - bound_pert,
                                        'center_y':self.y_pert - bound_pert}) 
        kwargs_upper_lens_light.append({'amp':0.,'R_sersic': 5, 
			                        	'n_sersic': 12,
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


        # Model parameters
        self.lens_params       = lens_params
        self.ps_params         = ps_params
        self.lens_light_params = lens_light_params

    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
	    from produce_kwargs_result import newProfile as produce_kwargs_result
	    return produce_kwargs_result(samples_mcmc,param_mcmc,i)
	    
    def str_param(self,param):  
        from get_param_title import newProfile as str_param
        return str_param(param)
