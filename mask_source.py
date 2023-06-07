import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip

from Utils.tools   import *
from Utils.get_res import *
from Data.input_data import init_kwrg_data
from Plots.plotting_tools import get_ModelPlot

def get_vlim(im):
    fact = 20.
    kw_vlim = {}
    dat = sigma_clip(im,sigma=3)
    vmax,vmin = np.max(dat),np.min(dat)
    kw_vlim["vmax"],kw_vlim["vmin"] = vmax*fact,vmin/fact
    return kw_vlim

def get_masksource(ref_setting,thr_source,thr_ps,save_mask=True,save_plot=False,stnd_masking=True):
    setting = get_setting_module(ref_setting,1)
    kwres   = get_kwres(setting_name=ref_setting)["kwargs_results"]
    
    modelplot = get_ModelPlot(setting=setting,kwargs_results=kwres)
    bandmodel = modelplot._band_plot_list[0]._bandmodel
    PointSources_image = bandmodel.image(**kwres,unconvolved=False,source_add=False,
                                      lens_light_add=False, point_source_add=True)
    Source_image       = bandmodel.image(**kwres,unconvolved=False,source_add=True,
                                      lens_light_add=False, point_source_add=False)
    
    mask = np.ones_like(Source_image)
    mask[np.where(Source_image>thr_source)]=0
    mask[np.where(PointSources_image>thr_ps)]=1
    
    if stnd_masking:
        _,stnd_mask = init_kwrg_data(setting,saveplots=False,return_mask=True)
        mask = mask*stnd_mask
    if save_mask:
        savename_mask = f"{get_savefigpath(setting)}/mask_source.json"
        save_json(data=mask,filename=savename_mask)
        print(f"Saved mask as {savename_mask}")
    if save_plot:
        kw_data = init_kwrg_data(setting,saveplots=False,return_mask=False)
        data = kw_data["image_data"]
        kw_vlim     = get_vlim(data)
        masked_data = data*mask
        fig,ax = plt.subplots(2,)
        ax[0].imshow(masked_data,origin="lower", cmap = plt.cm.gist_heat,**kw_vlim)
        ax[1].imshow(mask,origin="lower", cmap = plt.cm.gist_heat,**kw_vlim)
        plt.savefig(savename_mask.replace("json","png"))
        print(f"Saved mask as {savename_mask.replace('json','png')}")
    return mask


if __name__=="__main__":
    #############################
    present_program(sys.argv[0])
    #############################
    parser = argparse.ArgumentParser(description="The mask for the lensed source")
    parser.add_argument("-thr_s", "--threshold_source", type=int, dest="thr_s", default=.3,
                        help="Threshold for the pixel value of the source to be masked")
    parser.add_argument("-thr_ps", "--threshold_point_source", type=int, dest="thr_ps", default=4,
                        help="Threshold for the pixel value of the Point Source to NOT be masked")

    #parser.add_argument("-NP", "--no_plot",dest="no_plot", default=False,action="store_true",
    #                    help="Ignore the corner plots")
    parser.add_argument('SETTING_FILE',default=[],help="Setting file to which the mask is fitted")
    
    args = parser.parse_args()
    thr_source = args.thr_s
    thr_ps     = args.thr_ps
    setting_name =  args.SETTING_FILE
    get_masksource(ref_setting=setting_name,thr_ps=thr_ps,thr_source=thr_source,save_mask=True,save_plot=True)
    success(sys.argv[0])
