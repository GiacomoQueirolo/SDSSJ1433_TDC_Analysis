#!/usr/bin/env python
# coding: utf-8
# define a function mag that returns the magnification given the parameters


import sys
import argparse 

import corner
import numpy as np
from lenstronomy.LensModel.lens_model import LensModel

from Utils.tools import *
from Utils.order_images import image_order_mltf
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.get_res_multifilter import *
from Posterior_analysis.mag_remastered import labels,fr_nms,get_mag_ratio
from Data.Multiband_Data.input_data_mltf import init_lens_model_list_mltf
from Data.Multiband_Data.Param_mltf import conv_mcmc_i_to_kwargs_mltf,get_Param_mltf


@check_mltf_setting
def mag_mltf(multifilter_setting):
    multifilter_setting   = get_multifilter_setting_module(multifilter_setting)
    mcmc_mag_rt_file_path = multifilter_setting.get_savejson_path("mag_rt")
    mcmc_file_path        = multifilter_setting.get_savejson_path("mcmc_smpl")
                                                
    #First check if there is the mcmc_mag.json file; if so, load and return it, if not compute it and save it
    try:
        from os.path import getmtime as last_mod
        if last_mod(str(mcmc_mag_rt_file_path))<last_mod(str(mcmc_file_path)):
            print("Mag ratio file present, but outdated")
            raise
        mcmc_mag_ratio = load_whatever(mcmc_mag_rt_file_path)
        print("Mag ratio file already present")
        return mcmc_mag_ratio

    except:
        print(f"Warning: No {mcmc_mag_rt_file_path} found; We create one now.")
        pass
  
    ########
    # Load Samples_mcmc and Param_mcmc
    samples_mcmc = get_mcmc_smpl_mltf(multifilter_setting)
    mag_mcmc     = get_mag_mcmc_mltf(multifilter_setting,samples_mcmc)
    
    #I want to obtain the correct image order
    #########################################
    new_order = image_order_mltf(multifilter_setting,samples_mcmc)
    temp_mcmc = np.array(mag_mcmc).transpose()
    mcmc_i    = [temp_mcmc[i] for i in new_order] 
    mcmc_mag_ratio = (np.array(mcmc_i).transpose()).tolist() #shape = (n*steps, dimensions)
    
    multifilter_setting.savejson_data(mcmc_mag_ratio,mcmc_mag_rt_file_path)
     
    return mcmc_mag_ratio 


# MOD_LLFR


@check_mltf_setting
def get_mag_mcmc_mltf(multifilter_setting,samples_mcmc):
    lens_model_list  = init_lens_model_list_mltf(multifilter_setting)
    lens_model_class = LensModel(lens_model_list=lens_model_list)
    Param_class      = get_Param_mltf(multifilter_setting)
    mag_mcmc = []
    for i in range(len(samples_mcmc)):
        kwargs_result_i = conv_mcmc_i_to_kwargs_mltf(multifilter_sett=multifilter_setting,mcmc_i=samples_mcmc[i],Param_class=Param_class)
        kwargs_lens     = kwargs_result_i["kwargs_lens"]        
        kwargs_ps       = kwargs_result_i["kwargs_ps"]
        mag_ratio       = get_mag_ratio(lens_model_class,kwargs_lens,kwargs_ps)
        mag_mcmc.append(mag_ratio)
    return mag_mcmc




@check_mltf_setting
def flux_ratio_mltf(multifilter_setting,kwargs_result,outnames=False):
    # check in old/mag_remasterd_notes.py an old, correct while 
    # less precise and more complicated way to do it
    kwargs_ps = kwargs_result["kwargs_ps"][0]
    amp_i     = kwargs_ps["point_amp"]
    # Let's do it wrt image A
    #amp_max = np.max(amp_i)
    #FR = np.array(amp_i)/amp_max
    ####################
    # reorder them so that it is in alphabetical order
    samples_mcmc = get_mcmc_smpl_mltf(multifilter_setting)
    new_order    = image_order_mltf(multifilter_setting,samples_mcmc)
    amp_i        = [amp_i[i] for i in new_order] 
    ####################
    amp_A = amp_i[0]
    FR    = np.array(amp_i[1:])/amp_A
    if outnames:
        return FR,fr_nms
    else:
        return FR

 

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Obtain the magnification ratio mcmc chain and plot them - Multifilter version")
    parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-NFR","--NoFluxRatio", action="store_false", dest="FR", default=True,
                    help="Do NOT compute the expected Flux Ratio (wrt image A)")
    parser.add_argument('MULTIFILTER_SETTING',nargs="1",default=[],help="Multifilter_setting file to consider")
    args = parser.parse_args()
    
    multifilter_setting = args.MULTIFILTER_SETTING
    corner_plot         = args.corner_plot
    FR = args.FR
    if FR:
        from Data.input_data import *
        
    print("Analysing: ",multifilter_setting)
    
    mcmc_mag = mag_mltf(multifilter_setting,svpth=True)
    if corner_plot:
        plot = corner.corner(np.array(mcmc_mag), labels=labels, show_titles=True)
        plot.savefig(f"{multifilter_setting.savefig_path}/Mag.png")
    if FR:
        kwargs_result = get_kwres_mltf(multifilter_setting)["kwargs_results"]
        fr_i          = flux_ratio_mltf(multifilter_setting)        
        kw_fr         = {nm:fr for nm,fr  in zip(fr_nms,fr_i)}
        save_json(kw_fr,f"{multifilter_setting.savefig_path}/FR.json")
    success(sys.argv[0])

