#!/usr/bin/env python
# coding: utf-8

# # Modelling of the Quad SDSSJ144+6007 with HST image
# Copy from HST_HR.ipynb, modelling both _ws and not _ws setting file with the same program

import os,sys 
import pickle
import argparse
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

from lenstronomy.Plots.model_plot import ModelPlot 

from Utils import get_res
from Utils.tools import *
from Utils import rewrite_read_results
from Utils import last_commands
from Utils.check_success import check_success
from Posterior_analysis import mag_remastered
from Data.image_manipulation import *
#from custom_logL import logL_ellipticity_aligned as logL_ellipticity # MOD_CUSTOM_LIKE
#from Custom_Model.custom_logL import logL_ellipticity_qphi as  logL_ellipticity # MOD_CUSTOM_LIKE_II
from Data.input_data import init_kwrg_data,init_kwrg_psf,init_kwrg_numerics,get_kwargs_model
from Custom_Model.custom_logL import init_kwrg_custom_likelihood
from Posterior_analysis.source_pos import get_source_pos_MCMC


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(prog=sys.argv[0],description="Lens modelling program",formatter_class=CustomFormatter)
    parser.add_argument('-rt','--run_type',type=int,dest="run_type",default=0,help= """Type of run: 
        0 = standard, PSO_it = 800*rf   PSO_prt = 300*rf    MCMCb = 2000*rf  MCMCr = 8000*rf 
        1 = append  MCMCb = 2000*rf (only used if no previous MCMC found )  MCMCr=8000*rf 
        2 = test run  PSO_it = 3   PSO_prt = 3   MCMCb = 1 MCMCr = 2 
        3 = append test   MCMCb = 2  MCMCr=2 
    (PSO_it: PSO iterations, PSO_prt: PSO particles, MCMCr: MCMC run steps, MCMCb: MCMC burn in steps)\n""")
    parser.add_argument('-rf','--run_factor',type=float,dest="run_factor",default=2.,help="Run factor to have longer run")
    parser.add_argument('-tc','--threadCount',type=int,dest="threadCount",default=150,help="Number of CPU threads for the MCMC parallelisation (max=160)")
    parser.add_argument('SETTING_FILE',default="",help="setting file to model")

    args         = parser.parse_args()
    run_type     = args.run_type
    run_fact     = args.run_factor
    setting_name = get_setting_name(args.SETTING_FILE).replace(".py","")
    threadCount  = args.threadCount  
    RND = False #set a random start of the PSO
    n_run_cut = 50  # to re-implement
    #Model PSO/MCMC settings
    append_MC=False
    if run_type==0:
        n_iterations = int(700*run_fact) #number of iteration of the PSO run
        n_particles  = int(400*run_fact) #number of particles in PSO run
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

    np.seterr(all="ignore")


    backup_path   = "backup_results"
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    savefig_path  = get_savefigpath(setting_name,backup_path) 
    mkdir(savefig_path)
    mkdir(savemcmc_path)
    save_log_command(save_dir=savefig_path)
    backend_filename = get_backend_filename(setting_name)
    setting_path = find_setting_path(setting_name)
    os.system("cp "+setting_path+"/"+setting_name+".py "+savefig_path+".") #we copy the setting file to that directory
    sys.path.append(setting_path)


    setting = get_setting_module(setting_name).setting()
    CP = check_if_CP(setting)
    WS = check_if_WS(setting)

    if CP:
        print("WARNING: Considering the PEMD mass profile for the main lens")
    if WS:
        print("WARNING: this model DO NOT consider the source")

    if WS and setting_name[-3:]!="_ws":
        setting_dir=find_setting_path(setting_name)
        new_name = setting_dir+"/"+setting_name+"_ws.py"
        os.system("cp "+setting_dir+"/"+setting_name+".py "+new_name )
        line_prepender(new_name,"#SIMPLE COPY FROM "+setting_name+".py")
        setting_name=setting_name+"_ws"


    


    # datetime object containing current date and time
    now = datetime.now()
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    print("Running: ",sys.argv[0])
    print("Machine: ",os.uname()[1]) 
    print("Setting file:", setting_name)
    print("Started the :", dt_string)
 
    ##Printout of the results
    print_res = open(savefig_path+"results.txt","w")
    print_res.write("Results for "+setting_name+" \n")
    print_res.write("#################################\n\n")
    print_res.write("Date and time of start : "+ dt_string+"\n")
    print_res.write("Results obtained with "+sys.argv[0]+"\n")
    print_res.write("Setting file used settings/"+setting_name+"\n")
    print_res.write("Comments: "+setting.comments+"\n")
    print_res.write("append_MC: "+str(append_MC)+"\n")


    kwargs_data,mask = init_kwrg_data(setting,saveplots=True,backup_path=backup_path,return_mask=True)

    # ### PSF 
    kwargs_psf = init_kwrg_psf(setting,saveplots=True,backup_path=backup_path)


    # ## Modelling

    # Choice of profiles for the modelling
    # Lens (mass) profile with perturber and the Shear

    kwargs_params = {'lens_model':        setting.lens_params,
                    'point_source_model': setting.ps_params,
                    'lens_light_model':   setting.lens_light_params}
    if not WS:
        kwargs_params['source_model'] =   setting.source_params


    # ### Parameters for the PSO/MCMC runs
     
    kwargs_likelihood = init_kwrg_custom_likelihood(setting,mask,custom="qphi")
     

    ######
    kwargs_model            = get_kwargs_model(setting)
    lens_model_list         = kwargs_model['lens_model_list']
    lens_light_model_list   = kwargs_model['lens_light_model_list']
    point_source_model_list = kwargs_model['point_source_model_list']
    if not WS:
        source_model_list = kwargs_model['source_light_model_list']
    ######
        
    kwargs_numerics = init_kwrg_numerics(setting)
    multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    # if you have multiple  bands to be modeled simultaneously, you can append them to the mutli_band_list
    kwargs_data_joint = {'multi_band_list': multi_band_list, 
                         'multi_band_type': 'multi-linear'  
                         # 'multi-linear': every imaging band has independent solutions of the surface brightness, 
                         #'joint-linear': there is one joint solution of the linear coefficients \
                         # demanded across the bands.
                        }


    if setting.sub ==False:
        joint_lens_with_light=[[0,0,["center_x","center_y"]],[1,1,["center_x","center_y"]]]
    else:
        joint_lens_with_light=[[0,1,["center_x","center_y"]]]


    


    kwargs_constraints = {'num_point_source_list': [4], 
                          'solver_type': 'NONE',
                          'joint_lens_with_light':joint_lens_with_light}
    # mod free source
    if not WS and not setting.FS:
        kwargs_constraints['joint_source_with_point_source'] = [[0, 0]]
                         
    #  'joint_lens_with_light': list [[i_light, k_lens, ['param_name1', 'param_name2', ...]], [...], ...],
    #   joint parameter between lens model and lens light model
    mcmc_file_name = save_json_name(setting,savemcmc_path,filename="mcmc_smpl")
    mcmc_logL_file_name = save_json_name(setting,savemcmc_path,filename="mcmc_logL")
        
    if append_MC :
        try:
            try:
                mc_init_sample = np.array(get_res.load_whatever(mcmc_file_name))
            except:
                print("TEST: load hp5 backend file")
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
    from Custom_Model.my_lenstronomy.my_fitting_sequence import MyFittingSequence # ONLY IMPORT IT HERE OR IT BREAKS THE CODE

    fitting_seq = MyFittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints,\
                                  kwargs_likelihood, kwargs_params)
    if not append_MC:
        fitting_kwargs_list = [['MY_PSO', {'sigma_scale': 1., 'n_particles': n_particles, 
                                           'n_iterations': n_iterations,"path":savemcmc_path,"threadCount":threadCount}]]
    else:
        ####################################################################
        # if append True but no mcmc found, then we try to find the PSO and
        # start the MCMC from its best result. If no PSO, then we have to 
        # run it as usual
        ####################################################################
        fitting_kwargs_list =[]
        if mc_init_sample is None:
            from Custom_Model.init_mcmc_from_pso import create_mcmc_init
            mc_init_sample = create_mcmc_init(setting,backup_path=backup_path)
            if mc_init_sample is None:
                n_iterations = int(700*run_fact) #number of iteration of the PSO run PSO_it  = 700*rf
                n_particles  = int(400*run_fact) #number of particles in the PSO run PSO_prt = 400*rf 
                print("append is True, but no mcmc or PSO files were possible to be found/open. ")
                print(f"Starting standard run with PSO_it={n_iterations} , PSO_prt={n_particles}")

                fitting_kwargs_list = [['MY_PSO', {'sigma_scale': 1., 'n_particles': n_particles, 
                                               'n_iterations': n_iterations,"path":savemcmc_path,"threadCount":threadCount}]]

    if RND == False:
        np.random.seed(3)


    fitting_kwargs_list.append(['MCMC', {'n_burn': n_burn, 'n_run': n_run, 'walkerRatio': 10, 'sigma_scale': .1,\
                                         "threadCount":threadCount, 'init_samples':mc_init_sample,"backend_filename":backend_filename}])
    # First chain_list with only the PSO, the burn_in sequence and the first n_run_cut
    chain_list = fitting_seq.fit_sequence(fitting_kwargs_list) 
    sampler_type, mc_sample, param_mcmc, mc_logL  = chain_list[-1]
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
        save_mcmc_json(setting=setting,data=chain_list[0],filename="pso",backup_path=backup_path)
        
    kwargs_result   = fitting_seq.best_fit()
    param_file_name = savemcmc_path+setting_name.replace("settings","mcmc_prm")+".dat"
    with open(param_file_name, 'w+') as param_file:
        for i in range(len(param_mcmc)):
            param_file.writelines(param_mcmc[i]+",\n")


    # Reconstruct mcmc chain
    kw_mcmc        = get_res.get_mcmc(setting_name,backup_path)
    chain_list[-1] = ['MCMC',kw_mcmc["mcmc_smpl"],kw_mcmc["mcmc_prm"],kw_mcmc["mcmc_logL"]]
    samples_mcmc   = kw_mcmc["mcmc_smpl"]

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

    ##Printout of the results with errors
    get_res.get_sigma_kw(setting,mcmc_chain=chain_list[-1],print_res=print_res,save=True)

    # Plot of the obtained models
    v_min,v_max     = setting.v_min,setting.v_max
    res_min,res_max = setting.res_min,setting.res_max

    modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result,image_likelihood_mask_list=[mask.tolist()],\
                          arrow_size=0.02, cmap_string="gist_heat")

    if not WS:
        from Plots.plotting_tools import plot_model    as PM
    else:
        from Plots.plotting_tools import plot_model_WS as PM
        
    PM(modelPlot,savefig_path,v_min,v_max,res_min,res_max)
    
    #Printout of all results after obtaining the amplitude
    with open(savefig_path+"/read_results.data","wb") as f:
            pickle.dump(kwargs_result, f)
            
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


    #Normalised plot
    f, axes = plt.subplots(figsize=(10,7))
    modelPlot.normalized_residual_plot(ax=axes,v_min=res_min, v_max=res_max)
    plt.savefig(savefig_path+"normalised_residuals.png")
    plt.close()

    #Caustics
    f, axes = plt.subplots(figsize=(10,7))
    modelPlot.source_plot(ax=axes, deltaPix_source=0.01, numPix=1000, with_caustics=True)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig(savefig_path+"caustics.png")
    plt.close()

    f, axes = plt.subplots(figsize=(10,7))
    modelPlot.decomposition_plot(ax=axes, text='Point source position', source_add=False, \
                    lens_light_add=False, point_source_add=True, v_min=v_min, v_max=v_max)
    plt.savefig(savefig_path+"point_source_position.png")
    plt.close()

    #CHECK_FR

    #Since the time is the same for all images (considering no time delay, or negligible), we can consider the 
    # flux ratio to be amp_i/amp_max

    FR,ratio_name = mag_remastered.flux_ratio(setting,kwargs_result,outnames=True)
    print_res.write("Flux ratio for "+setting.filter_name+"\n")

    for i,FR_i in enumerate(FR):
        print_res.write("Flux ratio:"+str(FR_i)+" "+str(ratio_name[i])+"\n")    
    print_res.write("########################\n")


    # the results of the MCMC chain
    #MOD_SOURCE 
    kwargs_source,str_src = get_source_pos_MCMC(setting,svfg=True)
    print_res.write(str_src)

    #Closing result.txt
    print_res.close()

    last_kw = {"read_results":kwargs_result,
    #           "read_sigma_up":kwargs_sigma_upper,
    #           "read_sigma_low":kwargs_sigma_lower,
    #           "read_fermat":kwargs_fermat,
               "read_source":kwargs_source,
               "FR":FR}
    for nm in last_kw:
        pickle_results(last_kw[nm],nm,savefig_path=savefig_path)


    # add some final commands
        
    rewrite_read_results.rewrite_read_results(setting,cut_mcmc=0,backup_path=backup_path,save=True)

    for i in last_commands.progs:
        last_commands.last_command(setting_name, i,log=True,run=True) 

    check_success(setting_name,verbose=1)

    success(sys.argv[0])

