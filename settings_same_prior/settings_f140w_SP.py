#Initial data parameter settings for f140w with lens light subtracted
#####################################################################
import numpy as np
from conversion import conv_5 as conv
from same_prior import lens_prior
class setting():
    def __init__(self):
        self.setting_file_name = "f140w_SP"
        #Comment to insert in the backup file
        self.comments="Change wrt f140w_same_prior: change of lens prior and of structure of setting file, now as a class"   


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
        self.data_path="./data/F140W/"
        self.image_name="F140W_new_lens.fits"
        self.err_name="e.F140W_new_lens.fits"

        self.psf_name="small_sub_sym_psf_0001_F140W_qso_sub.fits"
        self.err_psf="e.small_sub_sym_psf_0001_F140W_qso_sub.fits"

        self.pssf = 1 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 15 #it should be the size of the non-supersampled kernel in pixel  #>note that this changed due to new psf
        self.lens_light_model_name ="F140W_llm.fits"
        #bool_bin = True #to bin the psf in the code
        self.bool_bin=False
        self.fixed_mag =[False] #Fixed magnification (one element per image)
        self.sub = True #subtracted light
        self.filter_name = "f140w"
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max=-1,1
        self.e_v_min,self.e_v_max =-1,0
        self.res_min, self.res_max= -4,4

        #Coordinates
        ############
        #CD matrix
        # from degree to arcsec
        CD1_1   = 1.04141997294039E-05*3600                                                  
        CD1_2   = 3.40688292375962E-05*3600                                                  
        CD2_1   = 3.40688292375962E-05*3600                                                  
        CD2_2   = -1.0414199729403E-05*3600

        self.transform_pix2angle = np.array([[CD1_1, CD1_2], [CD2_1, CD2_2]])  
        ra_at_xy_0,dec_at_xy_0 = 0,0

        # Parameter for the lens
        center_x_pix,center_y_pix = 39.129, 44.234
        #For the perturber
        x_pert_pix,y_pert_pix = 37, 29.4463
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [53.329, 35.7093, 25.3736, 41.718]
        y_im_pix = [33.4973, 32.8209, 42.252, 53.859]
        self.z_source = 2.737
        self.z_lens= 0.407 

        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image, self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale, self.ra_at_xy_0, self.dec_at_xy_0 = \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,ra_at_xy_0,dec_at_xy_0)

        #Mask coord
        ###########
        x_mask, y_mask =[65.94,66.62,76.22],[36.41,33.67,35.8] # coord in pixel in the HST image 
        #radius of the mask for the object(s) in pixels
        self.r_mask =  [2,3,2]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        #We implement the mask of the center
        x_mask_cen,y_mask_cen           =	np.array([40,41]),np.array([45,44])
        self.x_mask_cen,self.y_mask_cen =   x_mask_cen-1,y_mask_cen-1
        self.rad_mask_cen               =   [3,1]

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
			                     'center_y': self.center_y,
			                     'e1': lens_prior.e1,
		         	             'e2': lens_prior.e2})
        kwargs_lens_sigma.append({'center_x': .9*self.pix_scale, 
			                      'center_y': .9*self.pix_scale,
                                  'theta_E': lens_prior.sig_theta_E, 
                                  'e1': lens_prior.sig_e1, 
                                  'e2': lens_prior.sig_e2})

        kwargs_lower_lens.append({'center_x': self.center_x-8.*self.pix_scale,
                                  'center_y': self.center_y-8.*self.pix_scale,
                                  'theta_E': lens_prior.min_theta_E, 
                                  'e1': lens_prior.min_e1, 
                                  'e2': lens_prior.min_e2})

        kwargs_upper_lens.append({'center_x': self.center_x+8.*self.pix_scale, 
                                  'center_y': self.center_y+8.*self.pix_scale,
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
			                      'center_x': self.x_pert-2*self.pix_scale,
                                  'center_y': self.y_pert-2*self.pix_scale})

        kwargs_upper_lens.append({'theta_E': lens_prior.max_theta_pert, 
		                       	  'center_x': self.x_pert+2*self.pix_scale,
                                  'center_y': self.y_pert+2*self.pix_scale})

        #Shear
        fixed_lens.append({'ra_0':0., 'dec_0':0.}) 
        kwargs_lens_init.append({'gamma_ext' : lens_prior.gamma,     'psi_ext': lens_prior.psi})
        kwargs_lens_sigma.append({'gamma_ext': lens_prior.sig_gamma, 'psi_ext': lens_prior.sig_psi})
        kwargs_lower_lens.append({'gamma_ext': lens_prior.min_gamma, 'psi_ext': lens_prior.min_psi})
        kwargs_upper_lens.append({'gamma_ext': lens_prior.max_gamma, 'psi_ext': lens_prior.max_psi})
	        
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
        kwargs_lower_source = [{'amp': 0,'R_sersic': R_sersic_source - (sg_R_sersic_source*27) , 
				        'n_sersic': 0.5, 
				        'e1': e1_source -(sg_e1_source*27), 
				        'e2': e2_source - (sg_e2_source*27), 
				        'center_x': x_source - (sg_x_source*27), 
				        'center_y': y_source - (sg_y_source*27)}]

        kwargs_upper_source = [{'amp': 0,'R_sersic': R_sersic_source + (sg_R_sersic_source*27) , 
				        'n_sersic': 6., 
				        'e1': e1_source + (sg_e1_source*27), 
				        'e2': e2_source + (sg_e2_source*27), 
				        'center_x': x_source + (sg_x_source*27),
				        'center_y': y_source + (sg_y_source*27)}]

        #Ps light parameters 
        fixed_ps = [{}]
        lum_ps = np.array([0,0,0,0])
        kwargs_ps_init = [{'ra_image': self.x_image, 'dec_image': self.y_image,"point_amp": lum_ps}]
        bound_pix_ps=.5
        kwargs_ps_sigma = [{'ra_image' : 1*self.pix_scale*np.ones_like(self.x_image), 
                            'dec_image': 1*self.pix_scale*np.ones_like(self.x_image),
                            "point_amp":lum_ps}]
        kwargs_lower_ps = [{'ra_image' : self.x_image-bound_pix_ps*self.pix_scale,
                            'dec_image': self.y_image-bound_pix_ps*self.pix_scale,
                            "point_amp": lum_ps}]
        kwargs_upper_ps = [{'ra_image' : self.x_image+bound_pix_ps*self.pix_scale,
                            'dec_image': self.y_image+bound_pix_ps*self.pix_scale,
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
    					               'n_sersic':n_sersic_lens, 
					                   'center_x':self.x_pert,
					                   'center_y':self.y_pert})
        kwargs_lens_light_sigma.append({'amp':0.,'R_sersic': sg_R_sersic_lens*3,  
					                    'n_sersic': sg_n_sersic_lens*3, 
					                    'center_x':6.*self.pix_scale,
					                    'center_y':6.*self.pix_scale}) 
        kwargs_lower_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens - (sg_R_sersic_lens*9), 
					                    'n_sersic': 0.,
	                                    'center_x':self.x_pert - bound_pix_pert*self.pix_scale,
	                                    'center_y':self.y_pert - bound_pix_pert*self.pix_scale}) 
        kwargs_upper_lens_light.append({'amp':0.,'R_sersic':R_sersic_lens + (sg_R_sersic_lens*9), 
            					        'n_sersic':n_sersic_lens + (sg_n_sersic_lens*3),
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
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]


        # Model parameters
        self.lens_params = lens_params
        self.ps_params = ps_params
        self.lens_light_params = lens_light_params
        self.source_params = source_params


    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
        from produce_kwargs_result import standard as produce_kwargs_result
        return produce_kwargs_result(samples_mcmc,param_mcmc,i)
        

    def str_param(self,param):
        from get_param_title import standard as str_param
        return str_param(param)
