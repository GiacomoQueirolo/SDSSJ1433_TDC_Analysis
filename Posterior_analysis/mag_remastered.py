#!/usr/bin/env python
# coding: utf-8
# define a function mag that returns the magnification given the parameters
# In[1]:


import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
import os,sys
import pathlib as pth
import json
import argparse
import corner

from Utils.tools import *
from Utils.get_res import *
from Utils.order_images import get_new_image_order

labels = ["$\mu_B$/$\mu_A$","$\mu_C$/$\mu_A$","$\mu_D$/$\mu_A$"]
fr_nms = ["$FR_B/FR_A$","$FR_C/FR_A$","$FR_D/FR_A$"]
    
def mag(setting,svpth=False):
    setting = get_setting_module(setting).setting()
    #lnst_path     = pth.Path(setting_path).parent
    backup_path    = "./backup_results/"
    mcmc_dir_path  = get_savemcmcpath(setting,backup_path)
    mcmc_file_path = mcmc_dir_path+"/mcmc_smpl_"+strip_setting_name(setting)+".json"
    mcmc_mag_rt_file_path = pth.Path(str(mcmc_file_path).replace("smpl","mag_rt")) #mag ratio
    
    #this needed for the ordering of the images (if needed)
    savefig_path = get_savefigpath(setting,backup_path)
    
    
    #First check if there is the mcmc_mag.json file; if so, load and return it, if not compute it and save it
    try:
        from os.path import getmtime as last_mod
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

    
    
    ###########
    CP = check_if_CP(setting)
    ##########
    
    ########
    # Load Samples_mcmc and Param_mcmc
    kw_res       = get_mcmc(setting,backup_path)
    samples_mcmc = kw_res["mcmc_smpl"]
    param_mcmc   = kw_res["mcmc_prm"]
    ########
    
    mag_mcmc = get_mag_mcmc(samples_mcmc,param_mcmc,setting,CP)
    
    
    #I want to obtain the correct image order
    #########################################
    #read kwargs results
    kwargs_result = get_kwres(setting,backup_path)["kwargs_results"]

    # the first one is A no matter what
    
    #new_order = image_order(ra_im,dec_im)+1 #bc it gives the order assuming A in 0
    #new_order = [0,*new_order] 
    new_order = get_new_image_order(setting,samples_mcmc,starting_from_A=False)
    temp_mcmc = np.array(mag_mcmc).transpose()
    mcmc_i = [temp_mcmc[i] for i in new_order] 
    mcmc_mag_ratio = (np.array(mcmc_i).transpose()).tolist() #shape = (n*steps, dimensions)
    
    with mcmc_mag_rt_file_path.open(mode="w") as f:
        json.dump(mcmc_mag_ratio,f)
    if not svpth:
        return mcmc_mag_ratio
    else:
        return mcmc_mag_ratio, savefig_path


# In[ ]:


# MOD_LLFR
def get_mag_ratio(lens_model_class,kwargs_lens,kwargs_ps):
    cent_ra,cent_dec = kwargs_ps[0]["ra_image"] ,kwargs_ps[0]["dec_image"]
    mag = lens_model_class.magnification(*cent_ra,*cent_dec,kwargs_lens)
    # !! We assume correct order of the images !! #
    mag_ratio = mag[1:] / mag[0]
    return mag_ratio

def get_mag_mcmc(samples_mcmc,param_mcmc,setting,CP=False):
    if not CP:
        lens_model_list = ['SIE','SIS','SHEAR_GAMMA_PSI']
    else:
        lens_model_list = ['PEMD','SIS','SHEAR_GAMMA_PSI']
    lens_model_class = LensModel(lens_model_list=lens_model_list)
    mag_mcmc = []
    for i in range(len(samples_mcmc)):
        kwargs_result_i  = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
        kwargs_lens      = kwargs_result_i["kwargs_lens"]        
        kwargs_ps        = kwargs_result_i["kwargs_ps"]
        mag_ratio        = get_mag_ratio(lens_model_class,kwargs_lens,kwargs_ps)
        mag_mcmc.append(mag_ratio)
    return mag_mcmc


# In[ ]:


def flux_ratio(setting,kwargs_result,kwargs_numerics=None,kwargs_data=None,
               lens_model_list=None,light_model_list=None,point_source_list=['LENSED_POSITION'],outnames=False):
    setting=get_setting_module(setting,1)
    # check in old/mag_remasterd_notes.py an old, correct while 
    # less precise and more complicated way to do it
    kwargs_ps = kwargs_result["kwargs_ps"][0]
    amp_i     = kwargs_ps["point_amp"]
    # Let's do it wrt image A
    #amp_max = np.max(amp_i)
    #FR = np.array(amp_i)/amp_max
    ####################
    # reorder them so that it is in alphabetical order
    new_order    = get_new_image_order(setting,starting_from_A=True)
    amp_i        = [amp_i[i] for i in new_order] 
    ####################
    amp_A = amp_i[0]
    FR    = np.array(amp_i[1:])/amp_A
    if outnames:
        return FR,fr_nms
    else:
        return FR


# In[ ]:


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Obtain the magnification ratio mcmc chain and plot them")
    parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-NFR","--NoFluxRatio", action="store_false", dest="FR", default=True,
                    help="Do NOT compute the expected Flux Ratio (wrt image A)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    corner_plot = args.corner_plot
    FR = args.FR
    if FR:
        from input_data import *
        
    for sets in settings:
        if len(settings)>1:
            print("Analysing: ",sets)
        setting = get_setting_module(sets).setting()
        mcmc_mag,savefig_path = mag(sets,svpth=True)
        if corner_plot:
            plot = corner.corner(np.array(mcmc_mag), labels=labels, show_titles=True)
            plot.savefig(str(savefig_path)+"/Mag.png")
        if FR:
            kwargs_data     = init_kwrg_data(setting,saveplots=False)
            kwargs_numerics = init_kwrg_numerics(setting)
            kwargs_result   = get_kwres(sets)["kwargs_results"]
            
            #lens_model_list
            lens_model_list = ['SIE']
            if check_if_CP(setting):
                lens_model_list = ['PEMD']
            lens_model_list = [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']
            
            # light_model_list
            if setting.sub==False:
                light_model_list = ["SERSIC_ELLIPSE", "SERSIC","UNIFORM"]
            else:
                light_model_list = ["SERSIC","UNIFORM"]
            if hasattr(setting,"no_pert"):
                light_model_list=["UNIFORM"]
        
        
            fr_i  = flux_ratio(setting,kwargs_result,kwargs_numerics,kwargs_data,\
                             lens_model_list,light_model_list)
            
            kw_fr  = { nm:fr for nm,fr  in zip(fr_nms,fr_i)}
            with open(str(savefig_path)+"/FR.json","w") as f:
                json.dump(kw_fr,f)
    success(sys.argv[0])

