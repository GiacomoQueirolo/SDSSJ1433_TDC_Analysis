#Plot gif of image during n steps of the pso

import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import imageio

from tools import *
from get_res import *
from create_fits import *

from lenstronomy.Plots.model_plot import ModelPlot

from astropy.stats import sigma_clip
from Param import conv_mcmc_i_to_kwargs

def get_vlim(im,fact=20.):
    kw_vlim = {}
    dat = sigma_clip(im,sigma=3)
    vmax,vmin = np.max(dat),np.min(dat)
    kw_vlim["vmax"],kw_vlim["vmin"] = vmax*fact,vmin/fact
    return kw_vlim

def get_PSO_light(setting,n=10):
    setting     = get_setting_module(sett,1)
    pso_chain   = get_pso_chain(setting,backup_path=backup_path)
    lkl,pos,_   = pso_chain[1]
    bandmodel   = get_bandmodel(setting)
    kwdata,mask = init_kwrg_data(setting,return_mask=True)
    Orig        =  mask*kwdata["image_data"]
    
    kwargs_psf      = init_kwrg_psf(setting,saveplots=False)
    kwargs_numerics = init_kwrg_numerics(setting)
    multi_band_list = [[kwdata, kwargs_psf, kwargs_numerics]]
    kwargs_model    = get_kwargs_model(setting)
    kw_vlim         = get_vlim(Orig)

    ax = []
    for ind,pos_i in enumerate(pos) :
        if ind%n:
            # only consider one over n images
            continue
        kwres_i = conv_mcmc_i_to_kwargs(setting,pos_i)

        modelPlot_i = ModelPlot(multi_band_list, kwargs_model, kwres_i,image_likelihood_mask_list=[mask.tolist()],\
                      arrow_size=0.02, cmap_string="gist_heat")._band_plot_list[0]
    
        fig,axes = plt.subplots(1,4,figsize=(11,3))
        CL_Conv  = modelPlot_i._model #20*complete_model(setting,bandmodel=bandmodel,kwres=kwres_i,unconvolved=False)
        axes[0].text(1, 1, f"PSO It: {ind}", backgroundcolor="white",bbox=dict(fill=True, edgecolor='white', linewidth=2,color="white"))
        setting.v_min,setting.v_max= 0,3
        #print("setting.v_min",setting.v_min,"setting.v_max",setting.v_max)
        
        axes[0].imshow(CL_Conv,origin="lower", cmap = plt.cm.gist_heat,**kw_vlim )#vmin=setting.v_min,vmax=setting.v_max)
        axes[0].set_title("Model")
        axes[0].get_xaxis().set_visible(False)
        axes[0].get_yaxis().set_visible(False)
        axes[0].autoscale(False)
        axes[1].imshow(Orig,origin="lower", cmap = plt.cm.gist_heat,**kw_vlim) #vmin=setting.v_min,vmax=setting.v_max)
        axes[1].set_title("Data")
        axes[1].get_xaxis().set_visible(False)
        axes[1].get_yaxis().set_visible(False)
        axes[1].autoscale(False)

        Res      = modelPlot_i._norm_residuals#norm_residual(CL_Conv,kwdata["noise_map"],bandmodel)*mask      
        axes[3].imshow(Res,origin="lower", cmap = "bwr",vmin=-1.3,vmax=1.3) #vmin=setting.res_min,vmax=setting.res_max)

        axes[3].set_title("Residual")
        axes[3].get_xaxis().set_visible(False)
        axes[3].get_yaxis().set_visible(False)
        axes[3].autoscale(False)

        axes[2].set_title("-LogLikelihood")
        axes[2].plot(np.arange(len(lkl)),-np.array(lkl),color="g")
        axes[2].set_xlabel("Steps (PSO)")
        axes[2].set_ylabel("-LogLikelihood")
        iLklLim=int(len(lkl)/10)
        axes[2].set_ylim(min(-np.array(lkl)),-np.array(lkl)[iLklLim])
        axes[2].axvline(ind)
        axes[2].set_yscale('log')
        plt.tight_layout()
        fig.canvas.draw()
        data = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        data = data.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        ax.append(data)
        plt.close()
    return ax

if __name__=="__main__":
    ############################
    present_program(sys.argv[0]) 
    ############################
    parser = ArgumentParser(description="Create PSO gif of light")
    parser.add_argument("-fi","--fraction_of_index",default=10,dest="frac_ind", help="Which fraction of index(es) to consider (default=10)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting file(s) to model")
    args          = parser.parse_args()
    frac_ind      = int(args.frac_ind)
    setting_names = [s.replace(".py","") for s in get_setting_name(args.SETTING_FILES)]
    backup_path   = "backup_results"
    for sett in setting_names:
        setting     = get_setting_module(sett,1)
        print_setting(setting)
        axes = get_PSO_light(setting,n=frac_ind)
        seconds = 10
        imageio.mimsave(create_path_from_list([get_savefigpath(setting),"PSO_ConvLight.gif"]),axes,fps=len(axes)/seconds)
    success(sys.argv[0])
