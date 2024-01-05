#!/usr/bin/env python
# coding: utf-8

#### Def.source position from image position
import pickle
import argparse
import numpy as np
import pathos.multiprocessing as multiprocessing
from functools import partial
from corner import corner,quantile
from lenstronomy.LensModel.lens_model import LensModel

from Utils.tools import *
from Utils.get_res import *
from Data.input_data import init_lens_model

def pll_get_source_pos(kwargs_result_i,lensModel):
    ra_image   = kwargs_result_i["kwargs_ps"][0]["ra_image"]
    dec_image  = kwargs_result_i["kwargs_ps"][0]["dec_image"]
    source_x,source_y = get_source_gen(ra_image,dec_image,kwargs_result_i["kwargs_lens"],lensModel=lensModel)
    mc_i = []
    for j in range(len(source_x)):
        mc_i.append(source_x[j])
    for j in range(len(source_y)):
        mc_i.append(source_y[j])
    return mc_i

def test_get_source_pos_MCMC(setting,parallelised=True):
    #We implement the source position posterior
    #####################################
    setting      = get_setting_module(setting).setting()
    samples_mcmc = get_mcmc_smpl(setting)
    param_mcmc   = get_mcmc_prm(setting)
    labels_source = ["source position_x"]*4+["source position_y"]*4
    
    list_kwargs_results = []
    print(np.shape(samples_mcmc))
    samples_mcmc = samples_mcmc[:100]
    print(np.shape(samples_mcmc))
    for i in range(100):
        list_kwargs_results.append(setting.produce_kwargs_result(samples_mcmc,param_mcmc,i))
    lensModel = init_lens_model(setting)
    mcmc_source = []
    for i in range(len(samples_mcmc)):
        kwargs_result_i = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
        ra_image   = kwargs_result_i["kwargs_ps"][0]["ra_image"]
        dec_image  = kwargs_result_i["kwargs_ps"][0]["dec_image"]
        source_x,source_y = get_source_gen(ra_image,dec_image,kwargs_result_i["kwargs_lens"],setting)
        mc_i = []
        for j in range(len(source_x)):
            mc_i.append(source_x[j])
        for j in range(len(source_y)):
            mc_i.append(source_y[j])
        mcmc_source.append(mc_i)

    plot = corner(np.array(mcmc_source), labels=labels_source, show_titles=True)
    
    #then I want to consider the combined result:
    mcmc_source_T   = np.array(mcmc_source).transpose()
    mcmc_ra_source  = mcmc_source_T[0:4]
    mcmc_dec_source = mcmc_source_T[4:]

    mcmc_ra  = mcmc_ra_source.T.tolist()
    mcmc_dec = mcmc_dec_source.T.tolist()

    mcmc_combined_source = [np.hstack(np.transpose(mcmc_ra)),np.hstack(np.transpose(mcmc_dec))] 
    labels_source_comb = ["Combined source ra","Combined source dec"]
    plot = corner(np.transpose(mcmc_combined_source), labels=labels_source_comb, show_titles=True)

    # the number of points have to be reduced since now it's 4 times larger wrt the original sample size
    # due to the fact that we obtained them from 4 images
    resampled_combined = np.array([np.random.choice(i,len(samples_mcmc)) for i in mcmc_combined_source]).T
    #############################################
    # add the corner plot with all other params #
    #############################################
    print("np.shape(mcmc_combined_source),np.shape(samples_mcmc),np.shape(resampled_combined)")
    print(np.shape(mcmc_combined_source),np.shape(samples_mcmc),np.shape(resampled_combined))
    #mcmc_tot = np.array([*samples_mcmc,*resampled_combined ])
    mcmc_tot = np.hstack([samples_mcmc,resampled_combined])
    print("np.shape(mcmc_tot)") 
    print(np.shape(mcmc_tot)) 
    labels_tot = [*param_mcmc,*labels_source_comb]
    plot = corner(mcmc_tot, labels=labels_tot, show_titles=True)
    return 0    
        

def get_source_pos_MCMC(setting,svfg=False,output_mcmc=False,parallelised=True,samples_mcmc=None,param_mcmc=None,lensModel=None):
    #We implement the source position posterior
    #####################################
    setting      = get_setting_module(setting,1)
    if not samples_mcmc:
        samples_mcmc = get_mcmc_smpl(setting)
    if not param_mcmc:
        param_mcmc   = get_mcmc_prm(setting)
    labels_source = ["source position_x"]*4+["source position_y"]*4
    
    if parallelised:
        list_kwargs_results = []
        for i in range(len(samples_mcmc)):
            list_kwargs_results.append(setting.produce_kwargs_result(samples_mcmc,param_mcmc,i))
        if not lensModel:
            lensModel = init_lens_model(setting)
        with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
            mcmc_source = p.map(partial(pll_get_source_pos, lensModel=lensModel), list_kwargs_results)
    else:
        mcmc_source = []
        for i in range(len(samples_mcmc)):
            kwargs_result_i = setting.produce_kwargs_result(samples_mcmc,param_mcmc,i)
            ra_image   = kwargs_result_i["kwargs_ps"][0]["ra_image"]
            dec_image  = kwargs_result_i["kwargs_ps"][0]["dec_image"]
            source_x,source_y = get_source_gen(ra_image,dec_image,kwargs_result_i["kwargs_lens"],setting)
            mc_i = []
            for j in range(len(source_x)):
                mc_i.append(source_x[j])
            for j in range(len(source_y)):
                mc_i.append(source_y[j])
            mcmc_source.append(mc_i)
    if svfg:
        plot = corner(np.array(mcmc_source), labels=labels_source, show_titles=True)
        plot.savefig(get_savefigpath(setting)+"/MCMC_source.png")
    
    #then I want to consider the combined result:
    mcmc_source_T   = np.array(mcmc_source).transpose()
    mcmc_ra_source  = mcmc_source_T[0:4]
    mcmc_dec_source = mcmc_source_T[4:]

    mcmc_ra  = mcmc_ra_source.T.tolist()
    mcmc_dec = mcmc_dec_source.T.tolist()
        
    mcmc_combined_source = [np.hstack(np.transpose(mcmc_ra)),np.hstack(np.transpose(mcmc_dec))] 
    
    labels_source_comb = ["Combined source ra","Combined source dec"]
    if svfg:
        plot = corner(np.transpose(mcmc_combined_source), labels=labels_source_comb, show_titles=True)
        plot.savefig(get_savefigpath(setting)+"/MCMC_source_comb.png")

    #############################################
    # add the corner plot with all other params #
    #############################################
    if svfg:
        # the number of points have to be reduced since now it's 4 times larger wrt the original sample size
        # due to the fact that we obtained them from 4 images
        resampled_combined = np.array([np.random.choice(i,len(samples_mcmc)) for i in mcmc_combined_source]).T
        #mcmc_tot = np.transpose([*samples_mcmc, *mcmc_combined_source])
        mcmc_tot   = np.hstack([samples_mcmc,resampled_combined])
        labels_tot = [*param_mcmc,*labels_source_comb]
        plot = corner(mcmc_tot, labels=labels_tot, show_titles=True)
        plot.savefig(get_savefigpath(setting)+"/MCMC_post_ws.png")
    
        
    str_src = "\n#################################\n"
    str_src += "\nSource position\n" #Note that we are considering the position of the quasar too
    mcmc_source = [mcmc_ra,mcmc_dec]
    for i in range(2):
        val_min, val, val_max = quantile(np.array(mcmc_source[i]),q=[0.16, 0.5, 0.84])
        sig_min = np.abs(val_min-val)
        sig_max = val_max - val
        sig_ref = np.min([sig_max,sig_min])
        n_exp   = np.log10(1/sig_ref)
        fact    = pow(10,int(n_exp)+2)
        if sig_min==sig_max:
            str_src+=str(labels_source_comb[i]+"  " +str(np.trunc(np.array(val)*fact)/fact)+\
                            " +- "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
        else:
            str_src+=str(labels_source_comb[i]+" " +str(np.trunc(np.array(val)*fact)/fact)+\
                            " - "+str(np.trunc(np.array(sig_min)*fact)/fact)+\
                            " + "+str(np.trunc(np.array(sig_max)*fact)/fact)+"\n")
        if i==0:
            source_ra = (val,sig_min,sig_max)
        else:
            source_dec = (val,sig_min,sig_max)
    str_src+=str("\n#################################\n")        
    kwargs_source_out = {"source_ra":source_ra,"source_dec":source_dec}
    if output_mcmc:
        return kwargs_source_out,str_src,mcmc_source
    else:
        return kwargs_source_out,str_src

def get_source_pos_PSO(setting):
    setting = get_setting_module(setting,True)
    kwres = get_kwres(setting)["kwargs_results"]
    if not check_if_WS(setting):
        kw_source = kwres["kwargs_source"][0]
        center_x,center_y = kw_source["center_x"],kw_source["center_y"]
    else:
        ps_kw =  kwres["kwargs_ps"][0]
        ra_image,dec_image = ps_kw["ra_image"],ps_kw["dec_image"]
        print(ra_image,dec_image)
        source_x, source_y = get_source_gen(ra_image,dec_image,kwres["kwargs_lens"],setting)
        center_x, center_y = np.mean(source_x),np.mean(source_y)
    return {"source_ra": center_x,"source_dec":center_y}

def get_source_gen(ra_image,dec_image,kw_lens,setting=None,lensModel=None):
    if lensModel is None:
        lensModel = init_lens_model(setting)
    source_x, source_y = [],[]
    for j in range(len(ra_image)):
        x_source,y_source = lensModel.ray_shooting(ra_image[j],dec_image[j],kw_lens)
        source_x.append(x_source)
        source_y.append(y_source)
    return source_x,source_y

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Plot the posterior of the source position given the images positions from MCMC results (PSO also defined here)")
    parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-v","--verbose", action="store_true", dest="verbose", default=False,
                    help="Verbose")
    parser.add_argument("-pso","--PSO", action="store_true", dest="pso", default=False,
                    help="Print the PSO result too")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    corner_plot = args.corner_plot
    verbose = args.verbose
    pso = args.pso
    for sets in settings:
        if verbose:
            print(strip_setting_name(sets))
        kw_res, _ = get_source_pos_MCMC(sets,svfg=corner_plot)
        if corner_plot:
            print("Plot saved in "+get_savefigpath(sets))
        if verbose:
            print("MCMC result:\n",kw_res)
            if pso:
                print("PSO result:\n",get_source_pos_PSO(sets))
        with open(get_savefigpath(sets)+"/read_source.data","wb") as f:
            pickle.dump(kw_res, f)

    success(sys.argv[0])

