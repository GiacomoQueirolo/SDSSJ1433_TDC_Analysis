#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from tools import *
from get_res import *
from input_data import *
from image_manipulation import multifits_with_copied_hdr

import sys
from astropy.io import fits
from argparse import ArgumentParser
from lenstronomy.ImSim.MultiBand.single_band_multi_model import SingleBandMultiModel




def get_bandmodel(setting):
    kwargs_model     = get_kwargs_model(setting)
    kwargs_data,mask = init_kwrg_data(setting,return_mask=True)
    kwargs_numerics  = init_kwrg_numerics(setting)
    kwargs_psf       = init_kwrg_psf(setting,saveplots=False)
    multi_band_list  = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    image_likelihood_mask_list = [mask.tolist()]
    bandmodel = SingleBandMultiModel(multi_band_list, kwargs_model, likelihood_mask_list=image_likelihood_mask_list,band_index=0)
    return bandmodel

def source_light(setting,unconvolved=True,kwres=None,bandmodel=None):
    if kwres is None:
        kwres = get_kwres(setting)["kwargs_results"]
    if bandmodel is None:
        bandmodel = get_bandmodel(setting)
    lensed_source_model = bandmodel.image(unconvolved=unconvolved, 
        source_add=True,lens_light_add=False, point_source_add=False,**kwres)
    return lensed_source_model

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

if __name__=="__main__":
    present_program(sys.argv[0])
    parser = ArgumentParser(description="Create fits file with all light components")
    #parser.add_argument("-uc","--unconvolved",dest="unconvolved", default=False,action="store_true",
    #                help="Unconvolved source")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    settings    = get_setting_module(setting,sett=True)
    #unconvolved = args.unconvolved
    for sett in settings:
        image_file = sett.data_path+sett.image_name
        bandmodel  = get_bandmodel(sett)
        kwres      = get_kwres(sett)["kwargs_results"]
        SL  = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        LL  = lens_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        PSL = ps_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        CL  = complete_model(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=True)
        
        SL_Conv  = source_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        LL_Conv  = lens_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        PSL_Conv = ps_light(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        CL_Conv  = complete_model(sett,bandmodel=bandmodel,kwres=kwres,unconvolved=False)
        
        kwdata,mask = init_kwrg_data(sett,return_mask=True)
        Res         = mask*(CL_Conv - kwdata["image_data"])/kwdata["noise_map"]
        
        data_list  = [LL,SL,PSL,CL,LL_Conv,SL_Conv,PSL_Conv,CL_Conv,Res]
        data_title = ["Lens Light","Source Light","Point Source Light","Complete Light Model",\
        "Lens Light Convolved","Source Light Convolved","Point Source Light Convolved","Complete Light Model Convolved",\
                     "Residual"]
        
        fits_path = get_savefigpath(sett)+"/light_model.fits"

        multifits_with_copied_hdr(data_list=data_list,data_title=data_title,
                                  fits_parent_path=image_file,fits_res_namepath=fits_path,
                                  overwrite=True,verbose=True)
    success(sys.argv[0])

    

