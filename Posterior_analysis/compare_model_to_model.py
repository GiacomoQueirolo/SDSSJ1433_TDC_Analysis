# Compare_models
import sys
import argparse
import numpy as np
from copy import copy 
from functools import partial
from scipy.optimize import fmin

from Utils.tools import * 
from Utils.get_res import *
from Utils.create_fits import *
from Data.image_manipulation import create_mask
from Posterior_analysis.compare_model_to_data import create_model,to_min,get_bandmodel_given_model

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Fit a light model to another model. Takes 2 Setting. One is the result from which we take the model w.o. light amplitude, from the other we take the resulting model.\
                                     The amplitude of the model is then set to minimise the residual. This is then saved (only the absolute, the error is considered to be 1). ",
                        formatter_class=CustomFormatter,
                        usage='pokus --help',)
    parser.add_argument("-SM1", "--setting_model1", type=str,dest="sett_model1",  
                         help="Setting file from which the light model (apart the amplitude) is taken")
    parser.add_argument("-SM2", "--setting_model2", type=str,dest="sett_model2",  
                         help="Setting file from which the light model is taken to be compared to")
    ###############################
    present_program(sys.argv[0])
    ###############################
    args = parser.parse_args()
    # it's a copy from compare_model_to_data, we keep the same nomenclature: model1->model, model2->dat
    sett_model  = get_setting_module(args.sett_model1,1)
    model_str   = strip_setting_name(sett_model)
    sett_data   = get_setting_module(args.sett_model2,1)
    data_str    = strip_setting_name(sett_data,filter=True)
    backup_path = "./backup_results/" 
    res_dir     = f"{backup_path}/fitmodel2model/"
    res_dir     = f"{res_dir}/Mod_{model_str}_Data_{data_str}/"
    mkdir(res_dir)

    bandmod_res    = get_bandmodel_given_model(setting_data=sett_data,setting_model=sett_model)
    kwres_mod      = get_kwres(sett_model)["kwargs_results"]
    SL_bndmod_res  = source_light(sett_data,bandmodel=bandmod_res,unconvolved=False,kwres=kwres_mod)
    LL_bndmod_res  = lens_light(sett_data,bandmodel=bandmod_res,unconvolved=False,kwres=kwres_mod)
    PSL_bndmod_res = ps_light(sett_data,bandmodel=bandmod_res,unconvolved=False,kwres=kwres_mod)
    # The components which amplitude is free to vary are from the model, obtained as bandmodel
    # adapted to the grid of the "data" fitler
    kwargs_bndmod = {"SL":SL_bndmod_res,"LL":LL_bndmod_res,"PSL": PSL_bndmod_res}

    #kwmod_data     = init_kwrg_data(sett_data)
    #err_data       = kwmod_data["noise_map"]
    bandmod_data    = get_bandmodel(setting=sett_data)
    kwres_data      = get_kwres(sett_data)["kwargs_results"]
    SL_bndmod_data  = source_light(sett_data,bandmodel=bandmod_data,unconvolved=False,kwres=kwres_data)
    LL_bndmod_data  = lens_light(sett_data,bandmodel=bandmod_data,unconvolved=False,kwres=kwres_data)
    PSL_bndmod_data = ps_light(sett_data,bandmodel=bandmod_data,unconvolved=False,kwres=kwres_data)
    data            = SL_bndmod_data+LL_bndmod_data+PSL_bndmod_data
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

    data_pth  = f"{sett_data.data_path}/{sett_data.image_name}"
    err_pth   = f"{sett_data.data_path}/{sett_data.err_name}"

    fits_with_copied_hdr(model_min,fits_parent_path=data_pth,data_object=f"Light model of {model_str} fitted to light model {data_str}",
                        fits_res_namepath=f"{res_dir}/model.fits")
    
    fits_with_copied_hdr(np.abs(data-model_min)*mask_data,fits_parent_path=data_pth,data_object=f"Absolute Residual of Light model of {model_str} fitted to light model {data_str}",
                        fits_res_namepath=f"{res_dir}/abs_residual.fits")
    success(sys.argv[0])
