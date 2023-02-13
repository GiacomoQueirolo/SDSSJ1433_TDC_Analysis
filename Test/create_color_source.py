#!/usr/bin/env python
# coding: utf-8

# In[1]:


# from talk of 5th okt 22
# create a source image from the reconstructed one obtained from the different infrared filters
import sys
import numpy as np
from argparse import ArgumentParser

from Utils.tools import *
from Utils.get_res import *
from Data.input_data import *
from Utils.create_fits import source_light
from Data.conversion import conv_radec_to_xy
from Data.image_manipulation import fits_with_copied_hdr,shift_astrometry
from lenstronomy.ImSim.MultiBand.single_band_multi_model import SingleBandMultiModel


def get_bandmodel_2xSource(setting):
    kwargs_model     = get_kwargs_model(setting)
    kwargs_model["source_light_model_list"] = ["SERSIC","SERSIC"]
    kwargs_data,mask = init_kwrg_data(setting,return_mask=True)
    kwargs_numerics  = init_kwrg_numerics(setting)
    kwargs_psf       = init_kwrg_psf(setting,saveplots=False)
    multi_band_list  = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    image_likelihood_mask_list = [mask.tolist()]
    bandmodel = SingleBandMultiModel(multi_band_list, kwargs_model, likelihood_mask_list=image_likelihood_mask_list,band_index=0)
    return bandmodel

if __name__=="__main__":
    present_program(sys.argv[0])
    
    parser = ArgumentParser(description="Create fits file with delensed source in colors")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    settings    = get_setting_module(args.SETTING_FILES,sett=True)
    for i,sett in enumerate(settings):
        if not check_if_WS(sett):
            bandmodel = None
            if len(sett.source_params[0])==2:
                bandmodel = get_bandmodel_2xSource(sett)
            SL = source_light(sett,unconvolved=True,delensed=True,bandmodel=bandmodel)
            save_name = get_savefigpath(sett)+"/source_deconvlv_"+strip_setting_name(sett)+".fits"
            orig_imname = sett.data_path+sett.image_name
            fits = fits_with_copied_hdr(SL,orig_imname,data_object="Source deconvolved",fits_res_namepath=None)
            if i==0:
                A_radec   = [sett.x_image[0],sett.y_image[0]]
                coord_A   = np.array(A_radec) +np.array([218.3449834 , 60.12145745])
                coord_Axy = np.reshape(conv_radec_to_xy(sett,*A_radec),2)
                fits = shift_astrometry(fits,coord_A,coord_Axy)
            else:
                A_radec   = [sett.x_image[0],sett.y_image[0]]
                coord_Axy = np.reshape(conv_radec_to_xy(sett,*A_radec),2)
                fits = shift_astrometry(fits,coord_A,coord_Axy)
            print("Saving ",save_name)
            fits.writeto(save_name, overwrite=True)
    success(sys.argv[0])


# In[ ]:





# In[ ]:





# In[ ]:




