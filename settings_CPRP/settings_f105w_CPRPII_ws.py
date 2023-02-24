#SIMPLE COPY FROM settings_f105w_FINAL_wLL.py
#Initial data parameter settings for f105w with lens light subtracted
#####################################################################
import warnings
import numpy as np
from Data.conversion import conv_5 as conv 
from Data.conversion import find_barycenter
from Data.image_manipulation import get_transf_matrix
from same_prior_redone_I import lens_PEMD_prior as lens_prior
from same_prior_redone_I import *

class setting():

    warnings.simplefilter("ignore")

    def __init__(self):        
        #Comment to insert in the backup file
        self.comments="Change wrt settings_f105w_CPRPI_PSFTEST_ws: Error for the PSF corrected (a factor 5 too high)"
        #Modify here to indicate the image and the psf files and the model of lens light
        self.data_path   = "../old_gen_prog/data/F105W/"
        self.image_name  = "F105W_qso_skycorr.fits"
        self.err_name    = "e.F105W_qso_skycorr.fits"
        self.psf_name    = "test_psf_s5_F105W.fits"#"psf_SC_new_s5_F105W.fits"
        self.err_psf     = "e.test2_psf_s1_F105W.fits"#"e.psf_SC_new_s1_F105W.fits"
        self.filter_name = "f105w"
        self.pssf = 5 #point source supersampling factor -> how much the PSF was supersampled
        self.ssks = 31#it should be the size of the non-supersampled kernel in pixel
        #lens_light_model_name ="model_f814w_large.fits"
        self.lens_light_model_name = None
        self.bool_bin    = False
        self.fixed_mag   = [False] #Fixed magnification (one element per image)
        self.sub         = False #subtracted light
        self.already_sub = False
        #About this model
        self.WS = True #without source modelled
        self.CP = True #change Profile
        self.FS = False # MOD_Free_Source
        #"cosmetics": coutout for the printed images
        self.v_min,self.v_max     = -1,1
        self.e_v_min,self.e_v_max = -1,0
        self.res_min,self.res_max = -4,4

        #Coordinates
        ############
        self.transform_pix2angle          = get_transf_matrix(self.data_path+self.image_name,in_arcsec=True)
        self.ra_at_xy_0,self.dec_at_xy_0  = 0,0
    
        # Parameter for the lens
        center_x_pix,center_y_pix = 57.3934,57

        #For the perturber
        x_pert_pix, y_pert_pix = 54,41
        # Parameters for the source
        # order A (top), C (right), B (bottom),D (left)
        x_im_pix = [70.4763, 52.8799, 42.4159, 58.7773]
        y_im_pix = [45.4564, 44.7305, 54.0178, 65.5169]
        self.z_source = lens_prior.z_source
        self.z_lens   = lens_prior.z_lens

        # in ra,dec considering the rotation and deformation of the reference system
        self.x_image,self.y_image, self.x_pert, self.y_pert, self.center_x, self.center_y, self.pix_scale ,self.ra_at_xy_0, self.dec_at_xy_0= \
	        conv(x_im_pix,y_im_pix,center_x_pix,center_y_pix,x_pert_pix,y_pert_pix,self.transform_pix2angle,self.ra_at_xy_0,self.dec_at_xy_0)

        #Mask coord
        ###########
        x_mask, y_mask = [60,83,54,93,center_x_pix],[24,48,26,38,center_y_pix] # coord in pixel in the HST image 
        #radius of the mask for the object(s) in pixels
        self.r_mask = [3,3,1,2,2]
        self.x_mask,self.y_mask = np.array(x_mask)-1,np.array(y_mask)-1

        self.x_mask_cen,self.y_mask_cen = [],[]
        self.rad_mask_cen               = []
        self.R_tot_mask                 = 3./self.pix_scale # 3 "-> pixel
        self.new_center_mask_out        = find_barycenter(x_im_pix,y_im_pix)

        
        #Model parameters
        #################
        self.lens_params       = get_lens_params(self.center_x, self.center_y, self.x_pert, self.y_pert, self.pix_scale,CP=self.CP)
        self.ps_params         = get_ps_params(self.x_image,self.y_image,self.pix_scale)
        self.lens_light_params = get_lens_light_params(self.center_x,self.center_y,self.x_pert,self.y_pert,self.pix_scale,self.sub)
        
    #MOD_PROD
    def produce_kwargs_result(self,samples_mcmc,param_mcmc,i):
        from produce_kwargs_result import newProfile_with_main_lens_light as produce_kwargs_result
        return produce_kwargs_result(samples_mcmc,param_mcmc,i)
        
    def str_param(self,param):
        from get_param_title import newProfile_with_main_lens_light as str_param
        return str_param(param)
