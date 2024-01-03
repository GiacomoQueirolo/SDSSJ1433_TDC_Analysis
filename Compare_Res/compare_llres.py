# compare lens light profiles (p,phi and or e1,e2)

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

from lenstronomy.ImSim.MultiBand.single_band_multi_model import SingleBandMultiModel

from Utils.tools import *
from Utils.get_res import *
from Utils.create_fits import *
from Data.conversion import e12_from_qphi
from Data.input_data import get_kwargs_model
from Data.image_manipulation import fits_with_copied_hdr
from Utils.create_fits import lens_light,get_bandmodel

def diag_transfM(pixscale):
    return np.array([[pixscale,0],[0,pixscale]])

def roundup(f):
    if f-int(f)>0:
        return int(f)+1
    else:
        return int(f)


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    #parser = ArgumentParser(description="Plot comparison for ellipticity parameter of lens light results")
    parser = ArgumentParser(description="Plot comparison for ellipticity parameter of lens light results by plotting them")
    parser.add_argument("-n","--name_res",dest="name_res", help="Name of the result png",default="ellipt_comp.png")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting file(s) to model")
    args          = parser.parse_args()
    setting_names = [s.replace(".py","") for s in get_setting_name(args.SETTING_FILES)]
    backup_path   = "./backup_results"
    name_res      = backup_path+"/"+str(args.name_res).replace(".png","")+".png"
    
 
    f, axes = plt.subplots(3,roundup(len(setting_names)/3), figsize=(16, 16), sharex=False, sharey=False)
    kwargs_result = {"kwargs_lens_light":[{"center_x":1,"center_y":-2}]}
    dir_res = backup_path+"/tmp_llres/"
    mkdir(dir_res)
    for isett,sett in enumerate(setting_names):
        setting = get_setting_module(sett,1)
        print_setting(setting)
        bandmodel = get_bandmodel(setting=setting)
        
        if check_if_SUB(setting):
            q = setting.q_ll 
            phi = setting.phi_ll
            e1,e2 = e12_from_qphi(phi=phi,q=q)
            R_sersic,n_sersic = .75,6. # placeholders!!
            print("R and n for sersic profiles are set fix with placeholders for the optical images")
        else:
            kw_res_ll = get_kwres(setting)["kwargs_results"]["kwargs_lens_light"][0]
            R_sersic,n_sersic,e1,e2 = kw_res_ll["R_sersic"],kw_res_ll["n_sersic"],kw_res_ll["e1"],kw_res_ll["e2"]
        kwargs_model = get_kwargs_model(setting)
        kwargs_model["lens_light_model_list"] = ["SERSIC_ELLIPSE"]
        kwargs_model['lens_model_list'] =  None
        kwargs_model['point_source_model_list'] = None
        
        kwargs_result["kwargs_lens_light"][0]["R_sersic"] = R_sersic
        kwargs_result["kwargs_lens_light"][0]["n_sersic"] = n_sersic
        kwargs_result["kwargs_lens_light"][0]["e1"] = e1
        kwargs_result["kwargs_lens_light"][0]["e2"] = e2
        LL_Conv = lens_light(sett,bandmodel=bandmodel,kwres=kwargs_result,unconvolved=False)
        LL_nm   = strip_setting_name(setting)
        fits_path = create_path_from_list([".",setting.data_path,setting.image_name])
        fits_res_path = create_path_from_list([dir_res,LL_nm+".fits"])
        fits_res = fits_with_copied_hdr(data=LL_Conv,fits_parent_path=fits_path,data_object=LL_nm,fits_res_namepath=None)
        coord_A    = [218.3449834 , 60.12145745] # ra,dec
        coord_A_xy = np.array(conv_radec_to_xy(sett,kwargs_result["kwargs_lens_light"][0]["center_x"],kwargs_result["kwargs_lens_light"][0]["center_y"])).flatten()
        fits_corr  = shift_astrometry(fits_res,coord_A,coord_A_xy)
        fits_corr.writeto(fits_res_path, overwrite=True)
 
    success(sys.argv[0])
