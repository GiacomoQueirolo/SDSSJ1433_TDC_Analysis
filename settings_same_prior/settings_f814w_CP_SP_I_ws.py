#Initial data parameter settings for f814w with lens light subtracted
#####################################################################
import numpy  as np
from conversion import find_center
from conversion import conv_5 as conv
from same_prior_I import lens_PEMD_prior

class setting():
    def __init__(self):
    
        #Comment to insert in the backup file
        self.comments="Change wrt settings_f814w_CP_SP_ws.py:  use same_prior_I AND the new mask "
        #Modify here to indicate the image and the psf files and the model of lens light
        self.data_path="./data/F814W/"
        self.image_name="F814W_qso.fits"
        self.err_name="e.F814W_qso.fits"
        self.psf_name="psf_0005_F814W.fits"
        self.err_psf="e.psf_0001_F814W.fits"
        self.pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 41#it should be the size of the non-supersampled kernel in pixel
        self.lens_light_model_name=None
        self.bool_bin=False
        self.fixed_mag =[False] #Fixed magnification (one element per image)
        self.sub = True #subtracted light
        self.filter_name = "f814w"
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max=-4,1
        self.e_v_min,self.e_v_max =-3,0
        self.res_min, self.res_max= -3,3
        
        #Coordinates
        ############
        #CD matrix
        # from degree to arcsec

        CD1_1   = -8.1353962040008E-06*3600                                               
        CD1_2   = 7.41198940931543E-06*3600                                               
        CD2_1   = 7.41198940931543E-06*3600                                             
        CD2_2   = 8.13539620400086E-06*3600
         
        self.transform_pix2angle = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
        ra_at_xy_0,dec_at_xy_0 = 0,0

        # Parameter for the lens
        center_x_pix,center_y_pix = 62.6044,87.1731 

        #For the perturber
        x_pert_pix,y_pert_pix = 104.539,56
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [114.882,39.794,51.0978,92.899]
        y_im_pix = [109.799,103.96,39.6647,57.13]
        self.z_source = 2.737
        self.z_lens= 0.407
 
        #Conversion in angle in the cut image
        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0,self.dec_at_xy_0= \
	        conv (x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,ra_at_xy_0,dec_at_xy_0)

        #About this model
        self.WS = True #without source modelled
        self.CP = True #change Profile

        #Mask coord
        ###########
        x_mask, y_mask =[118],[148] # coord in pixel in the HST image 
        #radius of the mask for the object in pixels
        self.r_mask = [20.824865]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        #We implement the mask of the center (wrt the center, we need to add the numpix_cut/2  )
        x_mask_cen,y_mask_cen =	np.array([0,-1,+1,0,+1,0,-1])+center_x_pix,\
			        np.array([0,+1,-1,-2,-2,-1,0])+center_y_pix
        self.x_mask_cen,self.y_mask_cen = x_mask_cen-1,y_mask_cen-1
        self.rad_mask_cen               = [2,1,1,1,1,1,3]
        self.R_tot_mask                 = 3/self.pix_scale # 3 "-> pixel
        self.new_center_mask_out        = find_center(x_im_pix,y_im_pix)


        #Model parameters
        #################

        #Lens mass profile
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []

        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E': lens_PEMD_prior.theta_E,
                                 'center_x': self.center_x, 
			                     'center_y': self.center_y,
			                     'gamma': lens_PEMD_prior.gamma,
			                     'e1': lens_PEMD_prior.e1,
		         	             'e2': lens_PEMD_prior.e2})

        kwargs_lens_sigma.append({'center_x': .8*self.pix_scale, 
			                      'center_y': .8*self.pix_scale,
                                  'theta_E': lens_PEMD_prior.sig_theta_E, 
   			                      'gamma': lens_PEMD_prior.sig_gamma,
                                  'e1': lens_PEMD_prior.sig_e1, 
                                  'e2': lens_PEMD_prior.sig_e2})

        kwargs_lower_lens.append({'center_x': self.center_x-4.*self.pix_scale,
                                  'center_y': self.center_y-4.*self.pix_scale,
                                  'theta_E': lens_PEMD_prior.min_theta_E,
			                      'gamma': lens_PEMD_prior.min_gamma, 
                                  'e1': lens_PEMD_prior.min_e1, 
                                  'e2': lens_PEMD_prior.min_e2})

        kwargs_upper_lens.append({'center_x': self.center_x+6.*self.pix_scale, 
                                  'center_y': self.center_y+3.*self.pix_scale,
                                  'theta_E': lens_PEMD_prior.max_theta_E, 
			                      'gamma': lens_PEMD_prior.max_gamma,
                                  'e1': lens_PEMD_prior.max_e1, 
                                  'e2': lens_PEMD_prior.max_e2})

        # Perturber mass profile
        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E' : lens_PEMD_prior.theta_pert,
                                 'center_x': self.x_pert, 
                                 'center_y': self.y_pert})

        kwargs_lens_sigma.append({'theta_E': lens_PEMD_prior.sig_theta_pert, 
        			              'center_x': 2.*self.pix_scale,
                                  'center_y': 2.*self.pix_scale})

        kwargs_lower_lens.append({'theta_E': lens_PEMD_prior.min_theta_pert,
			                      'center_x': self.x_pert-6*self.pix_scale,
                                  'center_y': self.y_pert-6*self.pix_scale})

        kwargs_upper_lens.append({'theta_E': lens_PEMD_prior.max_theta_pert,
			                      'center_x': self.x_pert+6*self.pix_scale,
                                  'center_y': self.y_pert+6*self.pix_scale})

        #Shear
        fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
        kwargs_lens_init.append({'gamma_ext' : lens_PEMD_prior.gamma_ext,     'psi_ext': lens_PEMD_prior.psi_ext})
        kwargs_lens_sigma.append({'gamma_ext': lens_PEMD_prior.sig_gamma_ext, 'psi_ext': lens_PEMD_prior.sig_psi_ext})
        kwargs_lower_lens.append({'gamma_ext': lens_PEMD_prior.min_gamma_ext, 'psi_ext': lens_PEMD_prior.min_psi_ext})
        kwargs_upper_lens.append({'gamma_ext': lens_PEMD_prior.max_gamma_ext, 'psi_ext': lens_PEMD_prior.max_psi_ext})
        #Ps light parameters
        fixed_ps = [{}]
        #to consider the amplitude of the point source, we measure the maximum count in the image for each
        lum_ps = np.array([24.,22.,30.,5.])
        #Lower and upper boundary position freedom in pixel wrt the input value
        bound_pix_ps=.6 #1
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image': self.y_image,"point_amp": lum_ps}]
        kwargs_ps_sigma = [{'ra_image':  .2*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': .2*self.pix_scale*np.ones_like(self.x_image),"point_amp":lum_ps/10. }]
        kwargs_lower_ps = [{'ra_image':  self.x_image-bound_pix_ps*self.pix_scale, 
                            'dec_image': self.y_image-bound_pix_ps*self.pix_scale,"point_amp": .1*lum_ps}]
        kwargs_upper_ps = [{'ra_image':  self.x_image+bound_pix_ps*self.pix_scale, 
                            'dec_image': self.y_image+bound_pix_ps*self.pix_scale,"point_amp": 10*lum_ps}]

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
                                       'center_x':self.x_pert,'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':10.,'R_sersic':.1,  
				                        'n_sersic':.6, 
                                        'center_x':4.*self.pix_scale,'center_y':4.*self.pix_scale})
        kwargs_lower_lens_light.append({'amp':3,'R_sersic':0.1, 
				                        'n_sersic':1.,
                                        'center_x':self.x_pert-bound_pix_pert*self.pix_scale,
                                        'center_y':self.y_pert-bound_pix_pert*self.pix_scale})
        kwargs_upper_lens_light.append({'amp':25.,'R_sersic':1.,
			                            'n_sersic':6,
			                            'center_x':self.x_pert+bound_pix_pert*self.pix_scale,
			                            'center_y':self.y_pert+bound_pix_pert*self.pix_scale})

        #Uniform light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp':0.5})
        kwargs_lens_light_sigma.append({'amp':2})
        kwargs_lower_lens_light.append({'amp':0.})
        kwargs_upper_lens_light.append({'amp':5.})

        lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]


        # Model parameters
        self.lens_params = lens_params
        self.ps_params = ps_params
        self.lens_light_params = lens_light_params


    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
	    from produce_kwargs_result import newProfile
	    return newProfile(samples_mcmc,param_mcmc,i)
    def str_param(self,param):
        from get_param_title import newProfile
        return newProfile(param)
