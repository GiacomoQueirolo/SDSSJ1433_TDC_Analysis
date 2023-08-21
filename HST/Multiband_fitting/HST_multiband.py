#!/usr/bin/env python
# coding: utf-8

# Copy from HST_HRM_PLL.py
# adapted for multiple filters

import os,sys

import argparse

import numpy as np
import pickle
from datetime import datetime
import matplotlib.pyplot as plt

from lenstronomy.Data.psf import PSF
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Workflow.fitting_sequence import FittingSequence 

from Utils import get_res
from Utils.tools import *
from Utils.Multiband_Utils.Setting_multifilter import *
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.get_res_multifilter import get_mcmc_chain_mltf
#from Custom_Model.custom_logL import init_kwrg_custom_likelihood
from Custom_Model.Multiband_Model.custom_logL_mltf import init_kwrg_custom_likelihood_mltf
from Data.input_data import init_kwrg_data,init_kwrg_psf,init_kwrg_numerics 
from Data.Multiband_Data.input_data_mltf import get_kwargs_model_mltf,get_kwargs_constraints_mltf,get_fixed_sources_list_mltf
from Posterior_analysis.convergence_mcmc import test_convergence
from Plots.plotting_tools import plot_model 
from Plots.plotting_tools import plot_model_WS
import Posterior_analysis.Multiband_Posterior_Analysis.mag_mltf as mag_multifilter
from Posterior_analysis.Multiband_Posterior_Analysis.source_pos_mltf import get_source_pos_MCMC_mltf
from Utils.Multiband_Utils.rewrite_read_results_mltf import rewrite_read_results_mltf
"""
from Utils import last_commands 
from Utils.check_success import check_success
"""
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################

    parser = argparse.ArgumentParser(prog=sys.argv[0],description="Lens modelling program",formatter_class=CustomFormatter)
    parser.add_argument('-rt','--run_type',type=int,dest="run_type",default=0,help= """Type of run: 
        0 = standard, PSO_it = 1600*rf   PSO_prt = 600*rf    MCMCb = 1000*rf  MCMCr = 4000*rf 
        1 = append  MCMCb = 2000*rf (only used if no previous MCMC found )  MCMCr=8000*rf 
        2 = test run  PSO_it = 3   PSO_prt = 3   MCMCb = 1 MCMCr = 2 
        3 = append test   MCMCb = 2  MCMCr=2 
    (PSO_it: PSO iterations, PSO_prt: PSO particles, MCMCr: MCMC run steps, MCMCb: MCMC burn in steps)\n""")
    parser.add_argument('-rf','--run_factor',type=float,dest="run_factor",default=2.,help="Run factor to have longer run")
    parser.add_argument('-tc','--threadCount',type=int,dest="threadCount",default=150,help="Number of CPU threads for the MCMC parallelisation (max=160)")
    parser.add_argument('SETTING_FILES',nargs="+",default="",help="setting file to model")

    args        = parser.parse_args()
    run_type    = args.run_type
    run_fact    = args.run_factor
    settings    = get_setting_module(args.SETTING_FILES,1)
    threadCount = args.threadCount  
    RND = False #set a random start of the PSO
    n_run_cut = 50  # to re-implement
    #Model PSO/MCMC settings
    append_MC=False
    if run_type==0:
        n_iterations = int(800*run_fact) #number of iteration of the PSO run
        n_particles  = int(300*run_fact) #number of particles in PSO run
        n_run  = int(8000*run_fact) #MCMC total steps 
        n_burn = int(2000*run_fact) #MCMC burn in steps
    elif run_type ==1:
        append_MC=True
        # n_burn only used if previous MCMC not found
        n_burn = int(2000*run_fact) #MCMC burn in steps
        n_run  = int(8000*run_fact) #MCMC total steps 
    elif run_type==2:
        n_iterations = int(3) #number of iteration of the PSO run
        n_particles  = int(3) #number of particles in PSO run
        n_run  = int(2) #MCMC total steps 
        n_burn = int(1) #MCMC burn in steps
    elif run_type==3:
        append_MC   = True
        n_run  = int(2) #MCMC total steps 
        n_burn = int(1) #MCMC burn in steps
    else:
        raise RuntimeError("Give a valid run_type or implement it your own")

    np.seterr(all="ignore");


    backup_path      = "./backup_results_mltf"
    multifilter_sett = create_multifilter_setting(settings,backup_path=backup_path)
    mltf_sett_name   = get_multifilter_setting_name(multifilter_sett)
    savemcmc_path    = multifilter_sett.savemcmc_path
    savefig_path     = multifilter_sett.savefig_path
    save_log_command(save_dir=savefig_path)
    backend_filename = multifilter_sett.backend_filename
    for sett in settings:
        setting_name = get_setting_name(sett)    
        setting_path = find_setting_path(setting_name)
        os.system(f"cp {setting_path}/{setting_name} {savefig_path}/.") #we copy the setting file to that directory
    CP    = multifilter_sett.CP
    allWS = multifilter_sett.allWS
    if CP:
        print("WARNING: Considering the PEMD mass profile for the main lens")

    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print("Running: ",sys.argv[0])
    print("Machine: ",os.uname()[1]) 
    print("Multifilter_Setting file:",multifilter_sett.name )
    print("Started the :", dt_string)

    ##Printout of the results
    print_res = open(savefig_path+"results.txt","w")
    print_res.write("Results for "+mltf_sett_name+" \n")
    print_res.write("#################################\n\n")
    print_res.write("Date and time of start : "+ dt_string+"\n")
    print_res.write("Results obtained with "+sys.argv[0]+"\n")
    print_res.write("Multifilter Setting file used settings/"+mltf_sett_name+"\n")
    print_res.write("Comments: "+multifilter_sett.comments+"\n")
    print_res.write("append_MC: "+str(append_MC)+"\n")
    
    # Data
    multi_band_list = []
    masks           = []
    for sett in multifilter_sett.settings:
        # note: the plots are not saved
        kwargs_data,mask_i = init_kwrg_data(sett,saveplots=False,backup_path=backup_path,return_mask=True)
        ### PSF 
        kwargs_psf       = init_kwrg_psf(sett,saveplots=False,backup_path=backup_path)
        kwargs_numerics  = init_kwrg_numerics(sett)
        multi_band_list.append([kwargs_data, kwargs_psf, kwargs_numerics])
        masks.append(mask_i)
    kwargs_data_joint = {'multi_band_list': multi_band_list, 
                         'multi_band_type': 'multi-linear'  
                         # 'multi-linear': every imaging band has independent solutions of the surface brightness, 
                         #'joint-linear': there is one joint solution of the linear coefficients \
                         # demanded across the bands.
                        }

    # Modelling
    kwargs_model            = get_kwargs_model_mltf(multifilter_sett)
    lens_model_list         = kwargs_model['lens_model_list']
    lens_light_model_list   = kwargs_model['lens_light_model_list']
    point_source_model_list = kwargs_model['point_source_model_list']
    if not allWS:
        source_model_list = kwargs_model['source_light_model_list']    
    # Likelihood params
    kwargs_likelihood = init_kwrg_custom_likelihood_mltf(multifilter_sett,mask=masks,custom="PLL")
    # BE CAREFUL CHANGING THIS (MOD_PLL)
    ############################################################################################
    kwargs_constraints = get_kwargs_constraints_mltf(multifilter_sett,kwargs_model=kwargs_model)
    ############################################################################################

    # Params:
    # we assume that the frame of reference is good enough and thake the 1st setting lens_params as the lens params
    kwargs_params = {'lens_model':        multifilter_sett.lens_params,
                    'point_source_model': multifilter_sett.ps_params,
                    'lens_light_model':   multifilter_sett.lens_light_params}
    if not allWS:
        kwargs_params['source_model'] =  multifilter_sett.source_params
    ##########
            

    # "Hyper"Parameters for the PSO/MCMC runs
    mcmc_file_name      = multifilter_sett.get_savejson_path("mcmc_smpl")
    mcmc_logL_file_name = multifilter_sett.get_savejson_path("mcmc_logL")
    if append_MC :
        try:
            try:
                mc_init_sample = np.array(get_res.load_whatever(mcmc_file_name))
            except:
                import emcee
                mc_init_sample = np.array(emcee.backends.HDFBackend(backend_filename, read_only=True).get_chain(flat=True))
            n_burn = 0
        except:
            print(f"Both files {mcmc_file_name} and {backend_filename} not found or corrupted. Starting MCMC from scratch.")
            mc_init_sample = None
        try:
            mc_init_logL = get_res.load_whatever(mcmc_logL_file_name)
        except:
            mc_init_logL= np.array([])
    else:
        mc_init_sample = None
        mc_init_logL   = None

    # Try to solve the "OSError: [Errno 24] Too many open files" by deleting the 
    # n_run_cut implementation
    #from my_lenstronomy.my_fitting_sequence import MyFittingSequence # ONLY IMPORT IT HERE OR IT BREAKS THE CODE
    fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints,\
                                  kwargs_likelihood, kwargs_params)
    
    # in the first fitting we only use f475w as it's optical: high res, high S/N, low number of object to fit
    bands_compute_initial_step = [True if "475" in get_filter(sett) else False for sett in settings]
    other_bands = [True if "475" not in get_filter(sett) else False for sett in settings]
    all_bands = [True for _ in settings]
    update_settings_initial_step = {}
    update_settings_initial_step['kwargs_likelihood'] = {'bands_compute':bands_compute_initial_step}
    # fix all parameters that are only used in the other bands
    # lens and images will stay free
    # we can fix the sources as f475w has no source (not the positions)
    fixed_sources = get_fixed_sources_list_mltf(multifilter_sett)
    update_settings_initial_step['source_add_fixed'] = fixed_sources

    initial_fitting_kwlist = [['update_settings',update_settings_initial_step],
    ['PSO', {'sigma_scale': 1., 'n_particles': 200, 'n_iterations': 100,"threadCount":threadCount}],  # run PSO first band
    ['align_images', {'n_particles': 10*len(settings), 'n_iterations': 100*len(settings), 
                      'align_offset': True, 'align_rotation': True, 
                      'delta_shift': 0.2, 'delta_rot': 0.1,
                      'compute_bands': other_bands}],   
    ['update_settings', {'kwargs_likelihood': {'bands_compute': all_bands},  # now we fit all bands jointly
                        'source_remove_fixed': fixed_sources  # and release all parameters of the second band to be free
                        }]]

    if not append_MC:
        fitting_kwargs_list = [*initial_fitting_kwlist,['PSO', {'sigma_scale': 1., 'n_particles': n_particles, \
                            'n_iterations': n_iterations,"threadCount":threadCount}]]
    else:
        ####################################################################
        # if append True but no mcmc found, then we try to find the PSO and
        # start the MCMC from its best result. If no PSO, then we have to 
        # run it as usual
        ####################################################################
        fitting_kwargs_list = []
        if mc_init_sample is None:
            from Custom_Model.Multiband_Model.init_mcmc_from_pso_mltf import create_mcmc_init_mltf
            mc_init_sample = create_mcmc_init_mltf(multifilter_sett,backup_path=backup_path)
            if mc_init_sample is None:
                n_iterations = int(700*run_fact) #number of iteration of the PSO run PSO_it  = 700*rf
                n_particles  = int(400*run_fact) #number of particles in the PSO run PSO_prt = 400*rf 
                print("append is True, but no mcmc or PSO files were possible to be found/open. ")
                print(f"Starting standard run with PSO_it={n_iterations}, PSO_prt={n_particles}")

                fitting_kwargs_list = [*initial_fitting_kwlist,['PSO', {'sigma_scale': 1., 'n_particles': n_particles, 
                                    'n_iterations': n_iterations,"threadCount":threadCount}]]
        
    if RND == False:
        np.random.seed(3) 
    fitting_kwargs_list.append(['MCMC', {'n_burn': n_burn, 'n_run': n_run, 'walkerRatio': 10, 'sigma_scale': .1,\
                                         "threadCount":threadCount, 'init_samples':mc_init_sample,"backend_filename":backend_filename}])
    # First chain_list with only the PSO, the burn_in sequence and the first n_run_cut
    chain_list = fitting_seq.fit_sequence(fitting_kwargs_list) 
    sampler_type, mc_sample, param_mcmc, mc_logL  = chain_list[-1]
    # test mcmc convergence:
    test_convergence(setting=None,mcmc_logL=mc_logL)
    #append the previous results
    if append_MC:
        mc_sample = np.array([*mc_init_sample,*mc_sample])
        mc_logL   = np.array([*mc_init_logL,*mc_logL])

    # save here chain_list results
    save_json(data=mc_sample,filename=mcmc_file_name) 
    #save_mcmc_json(setting=setting,data=mc_sample, filename="mcmc_smpl",backup_path=backup_path)
    save_json(data=mc_logL,filename=mcmc_logL_file_name)
    #save_mcmc_json(setting=setting,data=mc_logL,   filename="mcmc_logL",backup_path=backup_path)
    if "PSO" in chain_list[0][0]:
        multifilter_sett.savejson_data(data=chain_list[0],filename="pso")
        
    kwargs_result   = fitting_seq.best_fit()
    param_file_name = savemcmc_path+multifilter_sett.get_savejson_name("prm").replace("json","dat")#multifilter_sett.savename.replace("mltf_setting","mcmc_prm").replace(".dll",".dat")
    with open(param_file_name, 'w+') as param_file:
        for i in range(len(param_mcmc)):
            param_file.writelines(param_mcmc[i]+",\n")

    # Reconstruct mcmc chain
    samples_mcmc,prm_mcmc,logL_mcmc =  [get_mcmc_chain_mltf(multifilter_sett,kw) for kw in ["smpl","prm","logL"]] #         = get_mcmc #multifilter_sett.get_mcmc()
    chain_list[-1] = ['MCMC',samples_mcmc,prm_mcmc,logL_mcmc]

    print_res.write("kwargs_model:"+str(kwargs_model)+"\n")
    print_res.write("kwargs_numerics:"+str(kwargs_numerics)+"\n")
    print_res.write("kwargs_constraints:"+str(kwargs_constraints)+"\n")
    del kwargs_likelihood["image_likelihood_mask_list"]
    print_res.write("kwargs_likelihood:"+str(kwargs_likelihood)+"\n")
    try:
        print_res.write("PSO particles: "+str(n_particles)+"\n")
        print_res.write("PSO run steps: "+str(n_iterations)+"\n")
    except NameError:
        pass 
    
    print_res.write("MCMC run steps: "+str(n_run)+"\n")
    print_res.write("MCMC burn in steps: "+str(n_burn)+"\n")
    print_res.write("number of non-linear parameters in the MCMC process: "+ str(len(param_mcmc))+"\n")
    print_res.write("parameters in order: "+str(param_mcmc)+"\n")
    print_res.write("number of evaluations in the MCMC process: "+str(np.shape(chain_list[-1][1])[0])+"\n")
    print_res.write("#################################\n")


    #Printout of all results after obtaining the amplitude
    with open(savefig_path+"/read_results.data","wb") as f:
            pickle.dump(kwargs_result, f)

    ##Printout of the results with errors -> for now ignored
    #get_res.get_sigma_kw(setting,mcmc_chain=chain_list[-1],print_res=print_res,save=True)

    # Plot of the obtained models
    multi_band_list_out = fitting_seq.multi_band_list
    modelPlot = ModelPlot(multi_band_list_out, kwargs_model, kwargs_result,image_likelihood_mask_list=masks,\
                          arrow_size=0.02, cmap_string="gist_heat")
    
    for band_index,sett in enumerate(settings):
        v_min,v_max     = sett.v_min,sett.v_max
        res_min,res_max = sett.res_min,sett.res_max
        if sett.WS:
            plot_model_WS(modelPlot,savefig_path,v_min,v_max,res_min,res_max,band=(band_index,get_filter(sett)))
        else:
            plot_model(modelPlot,savefig_path,v_min,v_max,res_min,res_max,band=(band_index,get_filter(sett)))
        #Normalised plot
        f, axes = plt.subplots(figsize=(10,7))
        modelPlot.normalized_residual_plot(ax=axes,band_index=band_index,v_min=res_min, v_max=res_max)
        plt.savefig(f"{savefig_path}/normalised_residuals_{get_filter(sett)}.png")
        plt.close()

    
    print_res.write("kwargs_results:\n")
    for res in kwargs_result:
        len_res = len(kwargs_result[str(res)])
        for i in range(len_res):
                print_res.write(str(res)+" "+str(i)+"\n")
                for j in kwargs_result[str(res)][i]:
                    print_res.write(str(j)+": "+str(np.trunc(np.array(kwargs_result[str(res)][i][str(j)])*1e4)/1e4)+"\n")
                print_res.write("\n")
        print_res.write("\n")

    print_res.write("\n#################################\n")
    logL   = modelPlot._imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
    n_data = modelPlot._imageModel.num_data_evaluate
    print_res.write(str(-logL * 2 / n_data)+' reduced X^2 of all evaluated imaging data combined\n')
    print_res.write("################################\n")

    #Caustics
    f, axes = plt.subplots(figsize=(10,7))
    modelPlot.source_plot(ax=axes, deltaPix_source=0.01, numPix=1000, with_caustics=True)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig(f"{savefig_path}/caustics.png")
    plt.close()

    f, axes = plt.subplots(figsize=(10,7))
    modelPlot.decomposition_plot(ax=axes, text='Point source position', source_add=False, \
                    lens_light_add=False, point_source_add=True)
    plt.savefig(f"{savefig_path}/point_source_position.png")
    plt.close()

    #CHECK_FR

    #Since the time is the same for all images (considering no time delay, or negligible), we can consider the 
    # flux ratio to be amp_i/amp_max

    FR,ratio_name = mag_multifilter.flux_ratio(multifilter_sett,kwargs_result,outnames=True)
    #print_res.write("Flux ratio for "+setting.filter_name+"\n")

    for i,FR_i in enumerate(FR):
        print_res.write("Flux ratio:"+str(FR_i)+" "+str(ratio_name[i])+"\n")    
    print_res.write("########################\n")


    # the results of the MCMC chain
    #MOD_SOURCE 
    kwargs_source,str_src = get_source_pos_MCMC_mltf(multifilter_sett,svfg=True)
    print_res.write(str_src)

    #Closing result.txt
    print_res.close()

    last_kw = {"read_results":kwargs_result,
               "read_source":kwargs_source,
               "FR":FR}
    for nm in last_kw:
        pickle_results(last_kw[nm],nm,savefig_path=savefig_path)

    rewrite_read_results_mltf(multifilter_sett,cut_mcmc=0,save=True)


    ## add some final commands
    #for i in last_commands.progs:
    #    last_commands.last_command(setting_name, i,log=True,run=True) 
    #check_success(setting_name,verbose=1)


    success(sys.argv[0])
