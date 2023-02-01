#!/usr/bin/env python
# coding: utf-8

# # Modelling of the Quad SDSSJ144+6007 with HST image
# Copy from HST_HR.ipynb, modelling both _ws and not _ws setting file with the same program
# In[1]:


import os,sys
import corner
import importlib
import numpy as np
import json,copy,pickle
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from lenstronomy.Data.psf import PSF
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource


# In[2]:


import get_res
from tools import *
import mag_remastered
from image_manipulation import *
#from custom_logL import logL_ellipticity_aligned as logL_ellipticity # MOD_CUSTOM_LIKE
from custom_logL import logL_ellipticity_qphi as logL_ellipticity # MOD_CUSTOM_LIKE_II
from input_data import init_lens_model,init_kwrg_data,init_kwrg_psf,init_kwrg_numerics,get_kwargs_model


# In[ ]:


############################
present_program(sys.argv[0])
############################


# In[ ]:


parser = ArgumentParser(description="Test of the PSO: run it many times and check how it goes")
parser.add_argument('-rt','--run_type',type=int,dest="run_type",default=0,help="Type of run: \n \
                0 = standard, PSO_it = 400*rf   PSO_prt = 200*rf   \
                (PSO_it: PSO iterations, PSO_prt: PSO particles, MCMCr: MCMC run steps, MCMCb: MCMC burn in steps)")
parser.add_argument('-rf','--run_factor',type=float,dest="run_factor",default=20.,help="Run factor to have longer run")
parser.add_argument('-nPSO','--num_PSO',type=int,dest="num_PSO",default=20,help="Number of PSOs: how many times we re-run the same PSO with different random seeds")
parser.add_argument('SETTING_FILE',default="",help="setting file to model")

args         = parser.parse_args()
run_type     = 0
run_fact     = args.run_factor
num_PSO    = args.num_PSO
setting_name = get_setting_name(args.SETTING_FILE).replace(".py","")
threadCount  = 160
RND = True #set a random start of the PSO
n_run_cut = 50  # to re-implement
#Model PSO/MCMC settings
append_MC=False
if run_type==0:
    n_iterations = int(400*run_fact) #number of iteration of the PSO run
    n_particles  = int(200*run_fact) #number of particles in PSO run
    n_run  = int(800*run_fact) #MCMC total steps 
    n_burn = int(200*run_fact) #MCMC burn in steps
elif run_type ==1:
    append_MC=True
    n_run  = int(800*run_fact) #MCMC total steps 
elif run_type==2:
    n_iterations = int(3) #number of iteration of the PSO run
    n_particles  = int(3) #number of particles in PSO run
    n_run  = int(2) #MCMC total steps 
    n_burn = int(1) #MCMC burn in steps
elif run_type==3:
    append_MC   = True
    n_run  = int(1) #MCMC total steps 
else:
    raise RuntimeError("Give a valid run_type or implement it your own")

np.seterr(all="ignore");


# In[3]:


backup_path   = "backup_results"
savemcmc_path = f"pso_test_{get_savemcmcpath(setting_name,backup_path)}"
savefig_path  = f"pso_test_{get_savefigpath(setting_name,backup_path)}"
setting_path = find_setting_path(setting_name)
mkdir(savefig_path)
mkdir(savemcmc_path)

setting = get_setting_module(setting_name).setting()

v_min,v_max     = setting.v_min,setting.v_max
res_min,res_max = setting.res_min,setting.res_max

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

# In[ ]:


kwargs_likelihood = {}
kwargs_likelihood["check_matched_source_position"] =True
kwargs_likelihood["source_position_tolerance"] = 0.01
kwargs_likelihood["force_no_add_image"] = True 
kwargs_likelihood["check_positive_flux"] = True  


# In[ ]:


# MOD_CUSTOM_LIKE
phi_ll = setting.phi_ll if setting.sub else None
q_ll   = setting.q_ll   if setting.sub else None
kwargs_likelihood["custom_logL_addition"] = logL_ellipticity(SUB=setting.sub,phi_ll=phi_ll,q_ll=q_ll)


# In[ ]:


kwargs_likelihood["image_likelihood_mask_list"] =  [mask.tolist()]


######
kwargs_model            = get_kwargs_model(setting)
lens_model_list         = kwargs_model['lens_model_list']
lens_light_model_list   = kwargs_model['lens_light_model_list']
point_source_model_list = kwargs_model['point_source_model_list']
if not WS:
    source_model_list  = kwargs_model['source_light_model_list']
######
    
kwargs_numerics  = init_kwrg_numerics(setting)
multi_band_list  = [[kwargs_data, kwargs_psf, kwargs_numerics]]
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

# Try to solve the "OSError: [Errno 24] Too many open files" by deleting the 
# n_run_cut implementation
from my_lenstronomy.my_fitting_sequence import MyFittingSequence  # ONLY IMPORT IT HERE OR IT BREAK THE CODE
from copy import copy

if not WS:
    from plotting_tools import plot_model    as PM
else:
    from plotting_tools import plot_model_WS as PM

def my_plot_mcmc_logL(ax, mcmc_logL, num_average=100):
        """
        plots the MCMC behaviour and looks for convergence of the chain
        :param samples_mcmc: parameters sampled 2d numpy array
        :param param_mcmc: list of parameters
        :param dist_mcmc: log likelihood of the chain
        :param num_average: number of samples to average (should coincide with the number of samples in the emcee process)
        :return:
        """
        num_samples = len(mcmc_logL)
        num_average = 1
        
        n_points = int((num_samples - num_samples % num_average) / num_average)
        mcmc_logL_averaged = -np.max(mcmc_logL[:int(n_points * num_average)].reshape(n_points, num_average), axis=1)
        mcmc_logL_normed   = (mcmc_logL_averaged - np.max(mcmc_logL_averaged)) / (np.max(mcmc_logL_averaged) - np.min(mcmc_logL_averaged))
        ax.scatter(np.arange(0,len(mcmc_logL_normed)),mcmc_logL_normed,color='k',marker=".", alpha=0.001)
        ax.set_title("Log Likelihood")
        ax.set_xlabel("MCMC Steps/"+str(num_average))
        return ax
def pickle_results(res,name,savefig_path):
    if not ".data" in name:
        name+=".data"
    with open(savefig_path+name,"wb") as f:
        pickle.dump(res, f)
    
for nPSO in range(num_PSO):
    np.random.seed(nPSO) # def random seed (all diff)
    savefig_path = f"{savefig_path}/pso_it_{nPSO}/"
    mkdir(savefig_path)
    
    os.system("cp "+setting_path+"/"+setting_name+".py "+savefig_path+".") #we copy the setting file to that directory

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


    fitting_seq = MyFittingSequence(copy(kwargs_data_joint), copy(kwargs_model), copy(kwargs_constraints),\
                                  copy(kwargs_likelihood), copy(kwargs_params))

    fitting_kwargs_list = [['MY_PSO', {'sigma_scale': 1., 'n_particles': n_particles, 
                                       'n_iterations': n_iterations,"path":savemcmc_path,"threadCount":threadCount}]]


    # First chain_list with only the PSO, the burn_in sequence and the first n_run_cut
    chain_list = fitting_seq.fit_sequence(fitting_kwargs_list) 


    # save here chain_list results
    pso_chain = chain_list[0]
    name_json = save_json_name(setting,savemcmc_path,"pso")
    save_json(pso_chain,name_json)
    
    kwargs_result   = fitting_seq.best_fit()
    
    f, axes = plt.subplots(1, 1, figsize=(18, 6))
    num_avrg = 100
    num_samples=len(pso_chain)
    if int(num_samples/num_avrg)>1000:
        num_avrg= num_samples/960
    my_plot_mcmc_logL(axes,pso_chain,num_avrg)
    del_ax = axes.plot([],[])
    del_ax[0].set_color("r")
    plt.savefig(savefig_path+'/MCMC_bhv_LogL.png')
    plt.close()

    print_res.write("kwargs_model:"+str(kwargs_model)+"\n")
    print_res.write("kwargs_numerics:"+str(kwargs_numerics)+"\n")
    print_res.write("kwargs_constraints:"+str(kwargs_constraints)+"\n")
    _kwargs_likelihood = copy(kwargs_likelihood)
    del _kwargs_likelihood["image_likelihood_mask_list"]
    print_res.write("kwargs_likelihood:"+str(_kwargs_likelihood)+"\n")
    print_res.write("PSO particles: "+str(n_particles)+"\n")
    print_res.write("PSO run steps: "+str(n_iterations)+"\n")
    print_res.write("#################################\n")

    modelPlot = ModelPlot(copy(multi_band_list), copy(kwargs_model), kwargs_result,likelihood_mask_list=[mask.tolist()],\
                          arrow_size=0.02, cmap_string="gist_heat")

    
    PM(modelPlot,savefig_path,v_min,v_max,res_min,res_max)

    pickle_results("read_results",kwargs_result,savefig_path)

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

    #Closing result.txt
    print_res.close()

success(sys.argv[0])

