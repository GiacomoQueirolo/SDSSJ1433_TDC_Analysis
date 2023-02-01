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

from util.tools_SDSS import *
from util.get_res_SDSS import *
from util.order_images_SDSS import get_new_image_order

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
    new_order = get_new_image_order(setting,starting_from_A=False)
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
        cent_ra,cent_dec = kwargs_ps[0]["ra_image"] ,kwargs_ps[0]["dec_image"]
        mag = lens_model_class.magnification(*cent_ra,*cent_dec,kwargs_lens)
        mag_ratio = mag[1:] / mag[0]
        mag_mcmc.append(mag_ratio)
    return mag_mcmc


# In[ ]:


def flux_ratio(setting,kwargs_result,kwargs_numerics,kwargs_data,
               lens_model_list,light_model_list,point_source_list=['LENSED_POSITION']):
    
    if setting.fixed_mag[0]:
        from lenstronomy.Data.imaging_data import ImageData
        from lenstronomy.Data.psf import PSF
        from lenstronomy.ImSim.image_model import ImageModel
        from lenstronomy.PointSource.point_source import PointSource
        from lenstronomy.LensModel.lens_model import LensModel
        from lenstronomy.LightModel.light_model import LightModel  
        
        lens_model_class = LensModel(lens_model_list=lens_model_list)

        #In order to do that we create the unconvolved modelled image of the ps alone, 
        # and consider the ratio of flux considering the "NONE" psf

        kwargs_lens = kwargs_result["kwargs_lens"]
        data_class = ImageData(**kwargs_data)
        psf_class  = PSF("NONE")
        kwargs_lens_light = kwargs_result["kwargs_lens_light"]
        lens_light_model_class = LightModel(light_model_list=light_model_list)

        kwargs_ps = kwargs_result["kwargs_ps"]
        # note: the relative magnification of point sources is not used as constraints in the fitting in the default settings of lenstronomy.
        # you can set this constraint with the keyword 'fixed_magnification_list' (see next block). 
        # The images are treated otherwise as separate linear amplitudes that are constraint independently of each other.
        point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[True])
        imageModel = ImageModel(data_class, psf_class, lens_model_class,
                                lens_light_model_class,point_source_class,
                                kwargs_numerics=kwargs_numerics)

        for i in range(len(kwargs_lens_light)):
            kwargs_lens_light[i]["amp"]=0

        image_ps = imageModel.image(kwargs_lens, kwargs_lens_light, kwargs_ps,unconvolved=False)

        
        cent_x = kwargs_result["kwargs_ps"][0]["ra_image"]/setting.pix_scale
        cent_y = kwargs_result["kwargs_ps"][0]["dec_image"]/setting.pix_scale

        rng_rd = 1#arcsec of the 1/2 square we consider
        rng    = rng_rd/setting.pix_scale
        numPix = kwargs_data["image_data"].shape[0] 
        amp_i  = []
        for j in range(len(cent_x)):
            x_i = cent_x[j] - rng +(numPix/2.) 
            x_f = cent_x[j] + rng +(numPix/2.) 
            y_i = cent_y[j] - rng +(numPix/2.) 
            y_f = cent_y[j] + rng +(numPix/2.) 
            im_i = []
            for i in range(len(image_ps)):
                im_line =[]
                if i>y_i and i<y_f:
                    for j in range(len(image_ps[i])):
                        if j>x_i and j<x_f:
                            im_line.append(image_ps[i][j])
                    im_i.append(im_line)
            im_i = np.array(im_i)
            amp_i.append(np.sum(im_i))
    else:
        kwargs_ps = kwargs_result["kwargs_ps"][0]
        amp_i=kwargs_ps["point_amp"]
    # Let's do it wrt image A
    #amp_max = np.max(amp_i)
    #FR = np.array(amp_i)/amp_max
    amp_A = amp_i[0]
    FR    = np.array(amp_i[1:])/amp_A
    return FR


# In[ ]:


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Obtain the magnification ratio mcmc chain and plot them")
    parser.add_argument("--corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-FR","--FluxRatio", action="store_false", dest="FR", default=True,
                    help="Compute the expected Flux Ratio (wrt image A)")
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

