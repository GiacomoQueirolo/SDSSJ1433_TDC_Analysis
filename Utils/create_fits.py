#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Utils.tools import *
from Utils.get_res import *
from Data.input_data import *
from Data.conversion import conv_radec_to_xy
from Data.image_manipulation import multifits_with_copied_hdr,shift_astrometry


import sys
import numpy as np
from astropy.io import fits
from argparse import ArgumentParser
from lenstronomy.ImSim.MultiBand.single_band_multi_model import SingleBandMultiModel

def get_bandmodel(setting):
    multi_band_list,mask = init_multi_band_list(setting=setting,return_mask=True,saveplots=False)
    image_likelihood_mask_list = [mask.tolist()]
    bandmodel = SingleBandMultiModel(multi_band_list, kwargs_model, likelihood_mask_list=image_likelihood_mask_list,band_index=0)
    return bandmodel

def source_light(setting,unconvolved=True,kwres=None,bandmodel=None,delensed=False):
    if kwres is None:
        kwres = get_kwres(setting)["kwargs_results"]
    if bandmodel is None:
        bandmodel = get_bandmodel(setting)
    if delensed:
        ra_grid, dec_grid = bandmodel.ImageNumerics.coordinates_evaluate
        source_model = bandmodel.SourceModel.surface_brightness(ra_grid, dec_grid, kwargs_source=kwres["kwargs_source"], k=None)
        source_model = bandmodel.ImageNumerics.re_size_convolve(source_model,unconvolved=unconvolved)
    else:
        source_model = bandmodel.image(unconvolved=unconvolved,source_add=True,lens_light_add=False, point_source_add=False,**kwres)
    return source_model

def lens_light(setting,unconvolved=True,kwres=None,bandmodel=None):
    if kwres is None:
        kwres = get_kwres(setting)["kwargs_results"]
    if bandmodel is None:
        bandmodel = get_bandmodel(setting)
    lens_light_model = bandmodel.image(unconvolved=unconvolved, 
        source_add=False,lens_light_add=True, point_source_add=False,**kwres)
    return lens_light_model

def ps_light(setting,unconvolved=True,kwres=None,bandmodel=None):
    if kwres is None:
        kwres = get_kwres(setting)["kwargs_results"]
    if bandmodel is None:
        bandmodel = get_bandmodel(setting)
    lensed_ps_model = bandmodel.image(unconvolved=unconvolved, 
        source_add=False,lens_light_add=False, point_source_add=True,**kwres)
    return lensed_ps_model

def complete_model(setting,unconvolved=True,kwres=None,bandmodel=None):
    if kwres is None:
        kwres = get_kwres(setting)["kwargs_results"]
    if bandmodel is None:
        bandmodel = get_bandmodel(setting)
    complete_model = bandmodel.image(unconvolved=unconvolved, 
         source_add=True,lens_light_add=True, point_source_add=True,**kwres)
    return complete_model

def norm_residual(model,error_map,bandmodel):
    norm_res = bandmodel.reduced_residuals(model, error_map=error_map)
    return norm_res

if __name__=="__main__":
    present_program(sys.argv[0])
    parser = ArgumentParser(description="Create fits file with all light components")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    settings = get_setting_module(args.SETTING_FILES,sett=True)
    for sett in settings:
        image_file = sett.data_path+sett.image_name
        bandmodel  = get_bandmodel(sett)
        kwres      = get_kwres(sett)["kwargs_results"]
        if not sett.WS:
            SL    = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
            SL_Dl = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True,delensed=True)
        LL  = lens_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        PSL = ps_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        CL  = complete_model(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        
        if not sett.WS:
            SL_Conv    = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
            SL_Conv_Dl = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False,delensed=True)
        LL_Conv  = lens_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        PSL_Conv = ps_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        CL_Conv  = complete_model(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        
        kwdata,mask  = init_kwrg_data(sett,return_mask=True)
        Orig         =  mask*kwdata["image_data"]
        #Res          =  mask*(CL_Conv - kwdata["image_data"])/kwdata["noise_map"]
        Res         = norm_residual(CL_Conv,kwdata["noise_map"],bandmodel)
        if not sett.WS:
            data_list   = [LL,SL,SL_Dl,PSL,CL,LL_Conv,SL_Conv,SL_Conv_Dl,PSL_Conv,CL_Conv,Orig,Res]
            data_object = ["Lens Light","Source Light","De-lensed Source Light","Point Source Light","Complete Light Model",\
            "Lens Light Convolved","Source Light Convolved","De-lensed Source Light Convolved","Point Source Light Convolved","Complete Light Model Convolved",\
            "Masked Original","Masked Normalised Residual"]
        else:
            data_list   = [LL,PSL,CL,LL_Conv,PSL_Conv,CL_Conv,Orig,Res]
            data_object = ["Lens Light","Point Source Light","Complete Light Model",\
            "Lens Light Convolved","Point Source Light Convolved","Complete Light Model Convolved",\
            "Masked Original","Masked Normalised Residual"]
        fits_path = get_savefigpath(sett)+"/light_model_"+sett.filter_name+".fits"

        multifits = multifits_with_copied_hdr(data_list=data_list,data_object=data_object,
                                  fits_parent_path=image_file,fits_res_namepath=None)
        # to be able to compare them, the astrometry has to be identical
        # unfortunately we are looking at slightly different astrometry
        # shifted by a certain ammount
        # I correct for that using image A as new reference point
        # its WCS corrdinates have been obtained from F814W and F475W (they have the right astrometry) using DS9
        # in agreement within ~3*10^-4 pixel or ~10^-5 arcec 
        coord_A    = [218.3449834 , 60.12145745] # ra,dec
        coord_A_xy = np.array(conv_radec_to_xy(sett,sett.ps_params[0][0]["ra_image"][0],sett.ps_params[0][0]["dec_image"][0])).flatten()
        mutlifits  = shift_astrometry(multifits,coord_A,coord_A_xy)
        mutlifits.writeto(fits_path, overwrite=True)
        print("Saving "+fits_path)
    success(sys.argv[0])

    

