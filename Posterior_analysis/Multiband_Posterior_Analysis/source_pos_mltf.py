#!/usr/bin/env python
# coding: utf-8


#### Def.source position from image position
import argparse
import numpy as np
from functools import partial
from corner import corner,quantile
import pathos.multiprocessing as multiprocessing


from Data.input_data import init_lens_model
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.get_res_multifilter import *
from Posterior_analysis.source_pos import pll_get_source_pos,get_source_gen
from Data.Multiband_Data.Param_mltf import conv_mcmc_i_to_kwargs_mltf,get_Param_mltf


"""
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
"""     

@check_mltf_setting
def get_source_pos_MCMC_mltf(multifilter_setting,svfg=False,output_mcmc=False,parallelised=True):
    #We implement the source position posterior
    #####################################
    
    samples_mcmc = get_mcmc_smpl_mltf(multifilter_setting)
    param_mcmc   = get_mcmc_prm_mltf(multifilter_setting)
    labels_source = ["source position_x"]*4+["source position_y"]*4
    Param_class   = get_Param_mltf(multifilter_setting)
    lensModel = init_lens_model(multifilter_setting.settings[0])
        
    if parallelised:
        list_kwargs_results = []
        for i in range(len(samples_mcmc)):
            kwres_i = conv_mcmc_i_to_kwargs_mltf(multifilter_sett=multifilter_setting,mcmc_i=samples_mcmc[i],Param_class=Param_class)
            list_kwargs_results.append(kwres_i)
        with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
            mcmc_source = p.map(partial(pll_get_source_pos, lensModel=lensModel), list_kwargs_results)
    else:
        mcmc_source = []
        for smpl_i in range(len(samples_mcmc)):
            kwargs_result_i = conv_mcmc_i_to_kwargs_mltf(multifilter_sett=multifilter_setting,mcmc_i=smpl_i,Param_class=Param_class)
            ra_image   = kwargs_result_i["kwargs_ps"][0]["ra_image"]
            dec_image  = kwargs_result_i["kwargs_ps"][0]["dec_image"]
            source_x,source_y = get_source_gen(ra_image,dec_image,kwargs_result_i["kwargs_lens"],lensModel=lensModel)
            mc_i = []
            for j in range(len(source_x)):
                mc_i.append(source_x[j])
            for j in range(len(source_y)):
                mc_i.append(source_y[j])
            mcmc_source.append(mc_i)
    if svfg:
        plot = corner(np.array(mcmc_source), labels=labels_source, show_titles=True)
        plot.savefig(f"{multifilter_setting.savefig_path}/MCMC_source.png")
    
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
        plot.savefig(f"{multifilter_setting.savefig_path}/MCMC_source_comb.png")

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
        plot.savefig(f"{multifilter_setting.savefig_path}/MCMC_post_ws.png")
    
        
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

@check_mltf_setting
def get_source_pos_PSO_mltf(multifilter_setting):
    kwres = get_kwres_mltf(multifilter_setting)["kwargs_results"]
    if not all([check_if_WS(sett) for sett in multifilter_setting.settings ]) :
        kw_source = kwres["kwargs_source"][0]
        center_x,center_y = kw_source["center_x"],kw_source["center_y"]
    else:
        ps_kw =  kwres["kwargs_ps"][0]
        ra_image,dec_image = ps_kw["ra_image"],ps_kw["dec_image"]
        print(ra_image,dec_image)
        source_x, source_y = get_source_gen(ra_image,dec_image,kwres["kwargs_lens"])
        center_x, center_y = np.mean(source_x),np.mean(source_y)
    return {"source_ra": center_x,"source_dec":center_y}

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Plot the posterior of the source position given the images positions from MCMC results (PSO also defined here) - MULTIFILTER VERSION")
    parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
                    help="DO NOT plot the corner plot")
    parser.add_argument("-v","--verbose", action="store_true", dest="verbose", default=False,
                    help="Verbose")
    parser.add_argument("-pso","--PSO", action="store_true", dest="pso", default=False,
                    help="Print the PSO result too")
    parser.add_argument('MULTIFILTER_SETTING',nargs="1",default=[],help="Multifilter_setting file to consider")
    args = parser.parse_args()
    
    multifilter_setting = get_multifilter_setting_module(args.MULTIFILTER_SETTING)
    corner_plot = args.corner_plot
    verbose = args.verbose
    pso = args.pso
    
    if verbose:
        print(strip_multifilter_name(multifilter_setting))
    kw_res, _ = get_source_pos_MCMC_mltf(multifilter_setting,svfg=corner_plot)
    if corner_plot:
        print(f"Plot saved in {multifilter_setting.savefig_path}")
    if verbose:
        print(f"MCMC result:\n{kw_res}")
        if pso:
            print(f"PSO result:\n{get_source_pos_PSO_mltf(multifilter_setting)}")
    success(sys.argv[0])

