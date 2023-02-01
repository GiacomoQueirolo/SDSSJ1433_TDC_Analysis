#Initial data parameter settings for f814w with lens light subtracted
#####################################################################
import warnings
import numpy as np
from same_prior import lens_prior
from util.conversion_SDSS  import conv_5, find_barycenter
from util.image_manipulation_SDSS import get_transf_matrix


class setting():

    warnings.simplefilter("ignore")

    def __init__(self):
        #Comment to insert in the backup file
        self.comments = "First setting file "
        #Modify here to indicate the image and the psf files and the model of lens light
        self.data_path   = "./data/"
        self.image_name  = "extr_r.fits"
        self.err_name    = "e.extr_r.fits"
        self.psf_name    = "psf_0005_r.fits"
        self.err_psf     = "e.psf_0001_r.fits"
        self.filter_name = "r"
        self.pssf = 5  #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 51 #it should be the size of the non-supersampled kernel in pixel
        self.bool_bin  = False
        self.fixed_mag = [False] #Fixed magnification (one element per image)
        
        self.sub = False #subtracted light
        self.lens_light_model_name = None   
        self.lens_light_model_name_full = None
        self.already_sub = False # is the lens light model already subtracted?
        
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max      = -4,1
        self.e_v_min,self.e_v_max  = -3,0
        self.res_min, self.res_max = -3,3

        #Coordinates
        ############
        #CD matrix
        self.transform_pix2angle         = get_transf_matrix(self.data_path+self.image_name,in_arcsec=True)
        self.ra_at_xy_0,self.dec_at_xy_0 = 0,0

        # Parameter for the lens
        center_x_pix,center_y_pix = 79.,46. 

        ##For the perturber
        #x_pert_pix,y_pert_pix = 104.538, 52
        x_forgr,y_foregr = 63.624,73.471
        # Parameters for the source
        # order A (top), B (bottom) , C (middle)
        x_im_pix = [84.727409, 69.707298, 79.]
        y_im_pix = [66.969235, 36.872343, 46.]

        self.z_source = 2.7887
        self.z_lens   = 0.76 #G4 , while G5 is 1.03 (http://arxiv.org/abs/1809.07337v1)
 
        #Conversion in angle in the cut image
        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image,self.x_forgr,self.y_forgr, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0,self.dec_at_xy_0= \
	        conv_5 (x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_forgr,y_foregr ,self.transform_pix2angle,self.ra_at_xy_0,self.dec_at_xy_0)
        
        
        # MOD_CUSTOM_LIKE_II
        #self.phi_ll = extract_phi_ll(model_name=self.lens_light_model_name_full,setting=self,min_ra=lens_prior.theta_E)#TMP find the correct model # extract_phi_ll(model_name="outdata_model_f814w.fits",setting=self,min_ra=lens_prior.theta_E)
        #self.phi_ll = extract_phi_ll(model_name=self.lens_light_model_name_full,setting=self,min_ra=lens_prior.theta_E)
        #self.q_ll   = extract_q_ll(model_name=self.lens_light_model_name_full,setting=self,min_ra=lens_prior.theta_E)

        #Mask coord
        ###########
        x_mask, y_mask = [36],[84] # coord in pixel in the HST image 
        #radius of the mask for the object in pixels
        self.r_mask = [10]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        #We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
        x_mask_cen,y_mask_cen           = [],[]
        self.x_mask_cen,self.y_mask_cen = [],[]#x_mask_cen-1,y_mask_cen-1
        self.rad_mask_cen               = []
        self.R_tot_mask                 = 32./self.pix_scale # 3 "-> pixel
        self.new_center_mask_out        = find_barycenter(x_im_pix,y_im_pix)

        #Model parameters
        #################

        #Lens mass profile
        bound_pix_mainlens = 10.
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
			                     'e1': lens_prior.e1,
		         	             'e2': lens_prior.e2})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_E, 
                                  'center_x': 5*self.pix_scale, 
			                      'center_y': 5*self.pix_scale,
                                  'e1': lens_prior.sig_e1, 
                                  'e2': lens_prior.sig_e2})

        kwargs_lower_lens.append({'theta_E': lens_prior.min_theta_E, 
                                  'center_x': self.center_x-bound_mainlens,
                                  'center_y': self.center_y-bound_mainlens,
                                  'e1': lens_prior.min_e1, 
                                  'e2': lens_prior.min_e2})

        kwargs_upper_lens.append({'theta_E': lens_prior.max_theta_E, 
			                      'center_x': self.center_x+bound_mainlens, 
                                  'center_y': self.center_y+bound_mainlens,
                                  'e1': lens_prior.max_e1, 
                                  'e2': lens_prior.max_e2})
                                  
        #Shear
        fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
        kwargs_lens_init.append({'gamma_ext' : lens_prior.gamma_ext,     'psi_ext': lens_prior.psi_ext})
        kwargs_lens_sigma.append({'gamma_ext': lens_prior.sig_gamma_ext, 'psi_ext': lens_prior.sig_psi_ext})
        kwargs_lower_lens.append({'gamma_ext': lens_prior.min_gamma_ext, 'psi_ext': lens_prior.min_psi_ext})
        kwargs_upper_lens.append({'gamma_ext': lens_prior.max_gamma_ext, 'psi_ext': lens_prior.max_psi_ext})
        
        #Ps light parameters
        fixed_ps = [{}]
        #Lower and upper boundary position freedom in pixel wrt the input value
        lum_ps   = np.array([0,0,0])
        bound_pix_ps = 2.
        bound_ps     = bound_pix_ps*self.pix_scale
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image':self.y_image,"point_amp": lum_ps}]
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),"point_amp": lum_ps}] #sigma pos was .6
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_ps, 'dec_image': self.y_image-bound_ps,"point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : self.x_image+bound_ps, 'dec_image': self.y_image+bound_ps,"point_amp": lum_ps}]
        
        #Lens light parameters
        fixed_lens_light        = []
        kwargs_lens_light_init  = []
        kwargs_lens_light_sigma = []
        kwargs_lower_lens_light = []
        kwargs_upper_lens_light = []
        #  light profile of foreground gal
        #Lower and upper boundary position freedom in pixel wrt the input value
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp':14.,'R_sersic':0.375, 
				                       'n_sersic':4.09, 
                                       'center_x':self.x_forgr,
                                       'center_y':self.y_forgr})
        kwargs_lens_light_sigma.append({'amp':10.,'R_sersic':.1,  
				                        'n_sersic':.6, 
                                        'center_x':3.*self.pix_scale,
                                        'center_y':3.*self.pix_scale})
        kwargs_lower_lens_light.append({'amp':3,'R_sersic':0.1, 
				                        'n_sersic':1.,
                                        'center_x':self.x_forgr-3.*self.pix_scale,
                                        'center_y':self.y_forgr-3.*self.pix_scale})
        kwargs_upper_lens_light.append({'amp':25.,'R_sersic':1.,
			                            'n_sersic':6,
			                            'center_x':self.x_forgr+3.*self.pix_scale,
			                            'center_y':self.y_forgr+3.*self.pix_scale})

        #Uniform light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp':0.5})
        kwargs_lens_light_sigma.append({'amp':2})
        kwargs_lower_lens_light.append({'amp':0.})
        kwargs_upper_lens_light.append({'amp':5.})

        lens_params       = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        ps_params         = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]


        # Model parameters
        self.lens_params       = lens_params
        self.ps_params         = ps_params
        self.lens_light_params = lens_light_params

    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
        from produce_kwargs_result import newProfile as produce_kwargs_result
        return produce_kwargs_result(samples_mcmc,param_mcmc,i)

    def str_param(self,param):
        from get_param_title import newProfile as str_param
        return str_param(param)
