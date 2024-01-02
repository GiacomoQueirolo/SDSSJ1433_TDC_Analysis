# Compare_models
import sys
import argparse
import numpy as np
from copy import copy 
from functools import partial
import matplotlib.pyplot as plt
from scipy.optimize import fmin

from Utils.tools import * 
from Utils.get_res import *
from Utils.create_fits import *
from Data.input_data import init_kwrg_data
from Data.image_manipulation import create_mask




def get_bandmodel_given_model(setting_data,setting_model):
    # we fit the setting_model model to the data of setting_data
    # e.g. the PSF is given by data
    kwargs_model     = get_kwargs_model(setting_model)
    kwargs_data,mask = init_kwrg_data(setting_data,return_mask=True)
    kwargs_numerics  = init_kwrg_numerics(setting_data)
    kwargs_psf       = init_kwrg_psf(setting_data,saveplots=False)
    multi_band_list  = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    image_likelihood_mask_list = [mask.tolist()]
    bandmodel = SingleBandMultiModel(multi_band_list, kwargs_model, likelihood_mask_list=image_likelihood_mask_list,band_index=0)
    return bandmodel

def create_model(SL,PSL,LL,fact_ps=1,fact_sl=1,fact_ll=1):
    return (SL*fact_sl)+(PSL*fact_ps)+(LL*fact_ll)

def diff_model(data,SL,PSL,LL,fact_ps=1,fact_sl=1,fact_ll=1):
    return np.abs(create_model(SL,PSL,LL,fact_ps,fact_sl,fact_ll)-data)

def to_min(x,mask_none,kwargs):
    fact_ps,fact_sl,fact_ll = x
    return np.nansum(mask_none*diff_model(fact_ps=fact_ps,fact_sl=fact_sl,fact_ll=fact_ll,**kwargs))

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Fit a light model to another data. Takes 2 Setting. One is the result from which we take the model w.o. light amplitude, from the other we take the observed light.\
                                     The amplitude of the model is then set to minimise the residual. This is then saved, both the absolute and the normalised one. ",
                        formatter_class=CustomFormatter,
                        usage='pokus --help',)
    parser.add_argument("-SM", "--setting_model", type=str,dest="sett_model",  
                         help="Setting file from which the light model (apart the amplitude) is taken")
    parser.add_argument("-SD", "--setting_data", type=str,dest="sett_data",  
                         help="Setting file from which the data is taken, to be compared to")
    ###############################
    present_program(sys.argv[0])
    ###############################
    args = parser.parse_args()
    
    sett_model  = get_setting_module(args.sett_model,1)
    model_str   = strip_setting_name(sett_model)
    sett_data   = get_setting_module(args.sett_data,1)
    data_str    = strip_setting_name(sett_data,filter=True)
    backup_path = "./backup_results/" 
    res_dir     = f"{backup_path}/fitmodel2data/"
    res_dir     = f"{res_dir}/Mod_{model_str}_Data_{data_str}/"
    mkdir(res_dir)

    bandmod       = get_bandmodel_given_model(setting_data=sett_data,setting_model=sett_model)
    kwres_mod     = get_kwres(sett_model)["kwargs_results"]
    SL_bndmod     = source_light(sett_data,bandmodel=bandmod,unconvolved=False,kwres=kwres_mod)
    LL_bndmod     = lens_light(sett_data,bandmodel=bandmod,unconvolved=False,kwres=kwres_mod)
    PSL_bndmod    = ps_light(sett_data,bandmodel=bandmod,unconvolved=False,kwres=kwres_mod)
    # The components which amplitude is free to vary are from the model, obtained as bandmodel
    # adapted to the grid of the "data" fitler
    kwargs_bndmod = {"SL":SL_bndmod,"LL":LL_bndmod,"PSL": PSL_bndmod}

    kwmod_data     = init_kwrg_data(sett_data)
    err_data       = kwmod_data["noise_map"]
    data           = kwmod_data["image_data"] 
    # the data to fit to is instead given by the setting_data
    kwargs_bndmod["data"] = data
    mask_data      = create_mask(data,setting=sett_data)
    mask_data_none = copy(mask_data)
    mask_data_none[np.where(mask_data==0)] = None 
    
    to_min_partial = partial(to_min,mask_none=mask_data_none,kwargs=kwargs_bndmod)
    # solve for the amplitude of the different component
    fct_ps_min,fct_sl_min,fct_ll_min = fmin(to_min_partial,[.5,1.,.7])
    # re-create the light model given the minimised amplitudes
    del kwargs_bndmod["data"]
    model_min = create_model(fact_ps=fct_ps_min,fact_sl=fct_sl_min,fact_ll=fct_ll_min,**kwargs_bndmod)
    # residual of the model given the error frame of the data
    residual  = norm_residual(model=model_min,error_map=err_data,bandmodel=bandmod)

    data_pth  = f"{sett_data.data_path}/{sett_data.image_name}"
    err_pth   = f"{sett_data.data_path}/{sett_data.err_name}"

    fits_with_copied_hdr(model_min,fits_parent_path=data_pth,data_object=f"Light model of {model_str} fitted to exposure {data_str}",
                        fits_res_namepath=f"{res_dir}/model.fits")

    fits_with_copied_hdr(residual,fits_parent_path=data_pth,data_object=f"Normalised Residual of Light model of {model_str} fitted to exposure {data_str}",
                        fits_res_namepath=f"{res_dir}/norm_residual.fits")

    fits_with_copied_hdr(np.abs(data-model_min)*mask_data,fits_parent_path=data_pth,data_object=f"Absolute Residual of Light model of {model_str} fitted to exposure {data_str}",
                        fits_res_namepath=f"{res_dir}/abs_residual.fits")
    success(sys.argv[0])
