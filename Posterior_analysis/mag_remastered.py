#!/usr/bin/env python
# coding: utf-8
# define a function mag that returns the magnification given the parameters

import corner
import argparse
import sys,json
import numpy as np
import multiprocess
import pathlib as pth
from os.path import getmtime as last_mod

from Utils.tools import *
from Utils.get_res import *
from Utils.order_images import get_new_image_order
from Posterior_analysis.tools_Post import get_array_BC
from Data.input_data import init_lens_model,init_kwrg_data,init_kwrg_numerics

labels_mr_BC = ["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_C$/$\mu_B$"]
labels_mr_AD = ["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_D$/$\mu_A$"]
labels_fr_BC = ["$FR_B/FR_A$","$FR_C/FR_A$","$FR_C/FR_B$"]
labels_fr_AD = ["$FR_B/FR_A$","$FR_C/FR_A$","$FR_D/FR_A$"]  
Warning_BC   = "Warning: we are considering the image couple BC instead of AD, this is not the standard analysis, be careful."

@check_setting
def gen_mag_ratio(setting,svpth=False,backup_path="./backup_results/",BC=False):
    mcmc_dir_path  = get_savemcmcpath(setting,backup_path)
    mcmc_file_path = save_json_name(setting,mcmc_dir_path,"mcmc_mag_rt") #mag ratio
    if BC:
        print(Warning_BC)
        mcmc_file_path = mcmc_file_path.replace("mcmc_mag_rt","mcmc_mag_rt_BC")
    mcmc_mag_rt_file_path = pth.Path(mcmc_file_path) 
    
    #this needed for the ordering of the images (if needed)
    savefig_path = get_savefigpath(setting,backup_path)
    
    #First check if there is the mcmc_mag.json file; if so, load and return it, if not compute it and save it
    try:
        if last_mod(str(mcmc_mag_rt_file_path))<last_mod(str(mcmc_file_path)):
            print("Mag ratio file present, but outdated")
            raise
        with mcmc_mag_rt_file_path.open() as f:
            mcmc_mag_ratio = json.load(f)
        print("Mag ratio file already present")
        if not svpth:
            return mcmc_mag_ratio
        else:
            return mcmc_mag_ratio,savefig_path
    except:
        print("Warning: No "+mcmc_mag_rt_file_path.name+" found; We create one now.")
        pass
    ########
    # Load Samples_mcmc and Param_mcmc
    kw_res       = get_mcmc(setting,backup_path)
    samples_mcmc = kw_res["mcmc_smpl"]
    param_mcmc   = kw_res["mcmc_prm"]
    ########
    #I want to obtain the correct image order
    #########################################
    # the first one is A no matter what
    new_order      = get_new_image_order(setting,samples_mcmc,starting_from_A=True,check_prev=True)
    mcmc_mag_ratio = gen_mag_ratio_mcmc(setting,samples_mcmc,param_mcmc,BC=BC,order=new_order)

    with mcmc_mag_rt_file_path.open(mode="w") as f:
        json.dump(mcmc_mag_ratio,f)
    if not svpth:
        return mcmc_mag_ratio
    else:
        return mcmc_mag_ratio, savefig_path
 
# MOD_LLFR
def gen_mag_ratio_i(lens_model_class,kwargs_lens,kwargs_ps,BC=False,order=None):
    mag = np.array(gen_mag_i(lens_model_class,kwargs_lens,kwargs_ps,order=order))
    mag_ratio = mag[1:] / mag[0]
    if BC:
        mag_ratio = get_array_BC(mag_ratio,relation="ratio")
    return np.array(mag_ratio).tolist()

def gen_mag_i(lens_model_class,kwargs_lens,kwargs_ps,order=None):
    cent_ra,cent_dec = kwargs_ps[0]["ra_image"] ,kwargs_ps[0]["dec_image"]
    mag = lens_model_class.magnification(cent_ra,cent_dec,kwargs_lens)
    if order is None:
        order = np.arange(len(mag))
    mag_ordered = [mag[i] for i in order] 
    return np.array(mag_ordered).tolist() 


@check_setting
def gen_mag_ratio_mcmc(setting,samples_mcmc=None,param_mcmc=None,order=None,BC=False,parall=True):
    if samples_mcmc is None:
        samples_mcmc = get_mcmc_smpl(setting_name=get_setting_name(setting))
    if param_mcmc is None:
        param_mcmc = get_mcmc_prm(setting_name=get_setting_name(setting))
    if order is None:
        order = get_new_image_order(setting=setting,mcmc=samples_mcmc,starting_from_A=True,check_prev=True)
    lens_model_class = init_lens_model(setting) 
    if parall:
        def _gen_mag_ratio_mcmc(i):
            kwargs_result_i  = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
            kwargs_lens      = kwargs_result_i["kwargs_lens"]        
            kwargs_ps        = kwargs_result_i["kwargs_ps"]
            mag_ratio        = gen_mag_ratio_i(lens_model_class,kwargs_lens,kwargs_ps,BC=BC,order=order)
            return mag_ratio
        with multiprocess.Pool() as pool:
            mag_mcmc  = pool.map(_gen_mag_ratio_mcmc,np.arange(0,len(samples_mcmc)))
    else:
        mag_mcmc = []    
        for i in range(len(samples_mcmc)):
            kwargs_result_i  = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
            kwargs_lens      = kwargs_result_i["kwargs_lens"]        
            kwargs_ps        = kwargs_result_i["kwargs_ps"]
            mag_ratio        = gen_mag_ratio_i(lens_model_class,kwargs_lens,kwargs_ps,order=order)
            mag_mcmc.append(mag_ratio)
    return mag_mcmc

@check_setting
def gen_mag_mcmc(setting,samples_mcmc=None,param_mcmc=None,order=None,parall=True):
    if samples_mcmc is None:
        samples_mcmc = get_mcmc_smpl(setting_name=get_setting_name(setting))
    if param_mcmc is None:
        param_mcmc = get_mcmc_prm(setting_name=get_setting_name(setting))
    if order is None:
        order = get_new_image_order(setting=setting,mcmc=samples_mcmc,starting_from_A=True,check_prev=True)
    lens_model_class = init_lens_model(setting) 
    if parall:
        def _gen_mag_mcmc(i):
            kwargs_result_i  = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
            kwargs_lens      = kwargs_result_i["kwargs_lens"]        
            kwargs_ps        = kwargs_result_i["kwargs_ps"]
            mag_ratio        = gen_mag_i(lens_model_class,kwargs_lens,kwargs_ps,order=order)
            return mag_ratio
        with multiprocess.Pool() as pool:
            mag_mcmc  = pool.map(_gen_mag_mcmc,np.arange(0,len(samples_mcmc)))
    else:
        mag_mcmc = []    
        for i in range(len(samples_mcmc)):
            kwargs_result_i  = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
            kwargs_lens      = kwargs_result_i["kwargs_lens"]        
            kwargs_ps        = kwargs_result_i["kwargs_ps"]
            mag_ratio        = gen_mag_i(lens_model_class,kwargs_lens,kwargs_ps,order=order)
            mag_mcmc.append(mag_ratio)
    return mag_mcmc

@check_setting
def flux_ratio(setting,kwargs_result,BC=False,outnames=False):
    # check in old/mag_remasterd_notes.py an old, correct while 
    # less precise and more complicated way to do it
    kwargs_ps = kwargs_result["kwargs_ps"][0]
    amp_i     = kwargs_ps["point_amp"]
    # Let's do it wrt image A 
    ####################
    # reorder them so that it is in alphabetical order
    new_order    = get_new_image_order(setting,starting_from_A=True)
    amp_i        = [amp_i[i] for i in new_order] 
    ####################
    amp_A = amp_i[0]
    FR    = np.array(amp_i[1:])/amp_A
    if BC:
        print(Warning_BC)
        FR = get_array_BC(FR,relation="ratio")

    if outnames:
        if BC:
            return FR,labels_fr_BC
        else:
            return FR,labels_fr_AD
    else:
        return FR


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Obtain the magnification ratio mcmc chain and plot them")
    parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-NFR","--NoFluxRatio", action="store_false", dest="FR", default=True,
                    help="Do NOT compute the expected Flux Ratio (wrt image A)")
    parser.add_argument("-BC", dest="BC", default=False,action="store_true",
                        help="Consider BC couple instead of AD (warning: not the standard)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings    = get_setting_module(args.SETTING_FILES,1)
    corner_plot = args.corner_plot
    FR          = args.FR
    BC          = args.BC 

    labels_mr = labels_mr_AD
    labels_fr = labels_fr_AD
    if BC:
        print(Warning_BC)
        labels_mr = labels_mr_BC
        labels_fr = labels_fr_BC

    for setting in settings:
        sets = get_setting_name(setting)
        if len(settings)>1:
            print("Analysing: ",sets)
        mcmc_mag,savefig_path = gen_mag_ratio(setting,svpth=True,BC=BC)
        if corner_plot:
            plot = corner.corner(np.array(mcmc_mag), labels=labels_mr, show_titles=True)
            plot.savefig(str(savefig_path)+"/Mag_Rt.png")
        if FR:
            kwargs_data     = init_kwrg_data(setting,saveplots=False)
            kwargs_numerics = init_kwrg_numerics(setting)
            kwargs_result   = get_kwres(sets)["kwargs_results"]     

            fr_i  = flux_ratio(setting,kwargs_result,BC=BC)
            
            kw_fr  = { nm:fr for nm,fr  in zip(labels_fr,fr_i)}
            with open(str(savefig_path)+"/FR.json","w") as f:
                json.dump(kw_fr,f)
    success(sys.argv[0])

