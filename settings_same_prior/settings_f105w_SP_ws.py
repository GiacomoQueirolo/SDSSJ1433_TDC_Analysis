#SIMPLE COPY FROM settings_f105w_FINAL_wLL.py
#Initial data parameter settings for f105w with lens light subtracted
#####################################################################
import numpy as np
from conversion import conv_5 as conv 
from same_prior import lens_prior
 
class setting():
    def __init__(self):        
        self.setting_file_name = "f105w_SP"
        #Comment to insert in the backup file
        self.comments="Change wrt to settings_f105w_same_prior: change of lens prior and of structure of setting file, now as a class"   

        R_sersic_lens = 0.11
        sg_R_sersic_lens = 0.05

        n_sersic_lens = 1.36 
        sg_n_sersic_lens= 0.5

        x_source = 0.466 
        y_source = -1.74
        sg_x_source = 0.02
        sg_y_source = 0.02


        #Modify here to indicate the image and the psf files and the model of lens light
        self.data_path="./data/F105W/"
        self.image_name="F105W_qso_lens_cut.fits"
        self.err_name="e.F105W_qso_lens_cut.fits"

        self.psf_name="big_sub_psf_0001_seA_F105W.fits"
        self.err_psf="e.big_sub_psf_0001_seA_F105W.fits"

        self.pssf = 1 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 17#it should be the size of the non-	supersampled kernel in pixel
        #lens_light_model_name ="model_f814w_large.fits"
        self.lens_light_model_name= None
        self.bool_bin=False
        self.fixed_mag =[False] #Fixed magnification (one element per image)
        self.sub = False #subtracted light
        self.filter_name = "f105w"
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max=-1,1
        self.e_v_min,self.e_v_max =-1,0
        self.res_min,self.res_max= -4,4

        #Coordinates
        ############
        #CD matrix
        # from degree to arcsec

        CD1_1   = 1.04151273131054E-05*3600  
        CD1_2   = 3.40685456789135E-05*3600                                                    
        CD2_1   = 3.40685456789135E-05*3600                                                    
        CD2_2   = -1.0415127313105E-05*3600  

        self.transform_pix2angle = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
        ra_at_xy_0, dec_at_xy_0 = 0,0
    
        # Parameter for the lens
        center_x_pix,center_y_pix = 57.3934,57

        #For the perturber
        x_pert_pix, y_pert_pix = 54,41
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [70.4763, 52.8799, 42.4159, 58.7773]
        y_im_pix = [45.4564, 44.7305, 54.0178, 65.5169]
        self.z_source = 2.737
        self.z_lens= 0.407


        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image,self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale ,self.ra_at_xy_0, self.dec_at_xy_0= \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,ra_at_xy_0,dec_at_xy_0)


        #Mask coord
        ###########
        x_mask, y_mask =[60,83,54,93,center_x_pix],[24,48,26,38,center_y_pix] # coord in pixel in the HST image 
        #radius of the mask for the object(s) in pixels
        self.r_mask = [3,3,1,2,2]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        self.x_mask_cen,self.y_mask_cen = [],[]
        self.rad_mask_cem               = []

        #Model parameters
        #################

        #Lens mass profile
        fixed_lens = []
        kwargs_lens_init = []
        kwargs_lens_sigma = []
        kwargs_lower_lens = []
        kwargs_upper_lens = []

        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E': lens_prior.theta_E,
                                 'center_x': self.center_x, 
			                     'center_y':self.center_y,
			                     'e1': lens_prior.e1,
		         	             'e2': lens_prior.e2})

        kwargs_lens_sigma.append({'center_x': .9*self.pix_scale, 
			                      'center_y': .9*self.pix_scale,
                                  'theta_E': lens_prior.sig_theta_E, 
                                  'e1': lens_prior.sig_e1, 
                                  'e2': lens_prior.sig_e2})

        kwargs_lower_lens.append({'center_x': self.center_x-4.*self.pix_scale,
                                  'center_y': self.center_y-2.*self.pix_scale,
                                  'theta_E': lens_prior.min_theta_E, 
                                  'e1': lens_prior.min_e1, 
                                  'e2': lens_prior.min_e2})

        kwargs_upper_lens.append({'center_x': self.center_x+2.*self.pix_scale, 
                                  'center_y': self.center_y+2.*self.pix_scale,
                                  'theta_E': lens_prior.max_theta_E, 
                                  'e1': lens_prior.max_e1, 
                                  'e2': lens_prior.max_e2})
                   
        # Perturber mass profile
        fixed_lens.append({})
        kwargs_lens_init.append({'theta_E' : lens_prior.theta_pert ,
                                 'center_x': self.x_pert, 
                                 'center_y': self.y_pert})

        kwargs_lens_sigma.append({'theta_E': lens_prior.sig_theta_pert , 
			                      'center_x': 6.*self.pix_scale,
                                  'center_y': 6.*self.pix_scale})

        kwargs_lower_lens.append({'theta_E': lens_prior.min_theta_pert,
			                      'center_x': self.x_pert-7*self.pix_scale,
                                  'center_y': self.y_pert-7*self.pix_scale})

        kwargs_upper_lens.append({'theta_E': lens_prior.max_theta_pert, 
		                       	  'center_x': self.x_pert+7*self.pix_scale,
                                  'center_y': self.y_pert+7*self.pix_scale})

        #Shear
        fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
        kwargs_lens_init.append({'gamma_ext' : lens_prior.gamma,     'psi_ext': lens_prior.psi})
        kwargs_lens_sigma.append({'gamma_ext': lens_prior.sig_gamma, 'psi_ext': lens_prior.sig_psi})
        kwargs_lower_lens.append({'gamma_ext': lens_prior.min_gamma, 'psi_ext': lens_prior.min_psi})
        kwargs_upper_lens.append({'gamma_ext': lens_prior.max_gamma, 'psi_ext': lens_prior.max_psi})
                          

        #Ps light parameters 
        fixed_ps = [{}]
        lum_ps = np.array([0,0,0,0])
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image': self.y_image,"point_amp": lum_ps}]
        bound_pix_ps=2.5
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),"point_amp":lum_ps}]
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_pix_ps*self.pix_scale, 
                            'dec_image': self.y_image-bound_pix_ps*self.pix_scale,"point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : self.x_image+bound_pix_ps*self.pix_scale, 
                            'dec_image': self.y_image+bound_pix_ps*self.pix_scale,"point_amp": lum_ps}]
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
		         		'center_x': self.center_x, 
		         		'center_y': self.center_y}]
        kwargs_lens_light_sigma=[{'amp': 0.,'R_sersic':.2, 
				        'n_sersic':0.5, 
				        'e1': .1,
				        'e2': 0.1,
				        'center_x': 0.8*self.pix_scale, 
				        'center_y': .8*self.pix_scale}]
        kwargs_lower_lens_light=[{'amp': 0.,'R_sersic': 1.,
				        'n_sersic':1.,
			         	'e1': .15 ,
			         	'e2': -.3,
				        'center_x': self.center_x-2.*self.pix_scale,
				        'center_y': self.center_y-4.*self.pix_scale,}]
        kwargs_upper_lens_light = [{'amp':0.,'R_sersic':5, 
				        'n_sersic':6,
			         	'e1': 0.5 ,
			         	'e2': 0.1,
				        'center_x': self.center_x+2.*self.pix_scale,
				        'center_y': self.center_y+2.*self.pix_scale}]
        # Perturber light profile
        bound_pix_pert=7 #5
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp' :0.,'R_sersic':R_sersic_lens, 
					                'n_sersic':n_sersic_lens , 
					                'center_x':self.x_pert,
					                'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					                    'n_sersic': sg_n_sersic_lens*3, 
					                    'center_x':6.*self.pix_scale,
					                    'center_y':6.*self.pix_scale}) 
        kwargs_lower_lens_light.append({'amp':0.,'R_sersic':.1, 
					                    'n_sersic': 1,
                                        'center_x':self.x_pert - bound_pix_pert*self.pix_scale,
                                        'center_y':self.y_pert - bound_pix_pert*self.pix_scale}) 
        kwargs_upper_lens_light.append({'amp':0.,'R_sersic':3, 
					                     'n_sersic': 8,
	                                     'center_x':self.x_pert + bound_pix_pert*self.pix_scale,
	                                     'center_y':self.y_pert + bound_pix_pert*self.pix_scale}) 

        #Uniform light profile
        fixed_lens_light.append({})
        kwargs_lens_light_init.append({'amp':1.5})
        kwargs_lens_light_sigma.append({'amp':3})
        kwargs_lower_lens_light.append({'amp':0.})
        kwargs_upper_lens_light.append({'amp':10.}) 
        
        lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        ps_params = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
        lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]
        

        # Model parameters
        self.lens_params = lens_params
        self.ps_params = ps_params
        self.lens_light_params = lens_light_params
        
    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
        from produce_kwargs_result import with_main_lens_light as produce_kwargs_result
        return produce_kwargs_result(samples_mcmc,param_mcmc,i)
        
    def str_param(self,param):
        from get_param_title import with_main_lens_light as str_param
        return str_param(param)
