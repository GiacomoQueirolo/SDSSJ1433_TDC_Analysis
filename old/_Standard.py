#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Mod from Notes 19th october: we consider a spline ml and see the difference with the only polynomial case 
#-> Expected: larger error distribution, smaller systematic


# In[2]:


"""
This script will find the generative noise model parameters that create mock lightcurves matching the data properties in term of gaussian and correlated noise
You can also provide directly the correct parameters in the config file. In this case I will just generate the python files to proceed to the step 3b and 3c
and skip the optimisation
""";


# In[ ]:


import copy
import os,sys
import glob,time
import importlib
import numpy as np
import pickle as pkl
import pathlib as pth
import argparse as ap
from corner import quantile
import matplotlib.pyplot as plt
from multiprocess import Pool, cpu_count
from matplotlib.ticker import FuncFormatter

#import pycs3.sim.run
import pycs3.gen.util
#import pycs3.sim.draw
import pycs3.gen.splml 
import pycs3.spl.topopt
import pycs3.gen.lc_func
import pycs3.sim.twk as twk
import pycs3.pipe.optimiser
import pycs3.pipe.pipe_utils as ut

####################
from tools import *
from stnd_plot import delayplot,dmagplot
from stnd_handling_data import * #Error,getresults,combine_series, combine_series_methodB
from plot_distrib import plot_err
from stnd_red_chi import get_chi_red
from pycs3_mod.sim.run import multirun
from pycs3_mod.sim.draw import multidraw
from My2Script import analyse_lc,analyse_lc_ij
from inspect_results import plt_err, plt_err_tot, plt_intr_err
####################


# In[4]:


#################################
########### Script 1 ############
#################################

def setup_directories(config):
    if not os.path.exists(config.config_directory):
        print("I will create the config directory for you ! ")
        mkdir(config.config_directory)
    if not os.path.exists(config.general_directory):
        print("I will create the general directory for you ! ")
        mkdir(config.general_directory)
    if not os.path.exists(config.lens_directory):
        print("I will create the lens directory for you ! ")
        mkdir(config.lens_directory)
    if not os.path.exists(config.analysis_directory):
        print("I will create the analysis directory for you ! ")
        mkdir(config.analysis_directory)   
    if not os.path.exists(config.simu_directory):
        print("I will create the simulation directory for you ! ")
        mkdir(config.simu_directory)
        
    if not os.path.exists(config.figure_directory):
        print("I will create the figure directory for you ! ")
        mkdir(config.figure_directory)
    if not os.path.exists(config.report_directory):
        print("I will create the report directory for you ! ")
        mkdir(config.report_directory)
    


# In[ ]:


#################################
########### Script 2 ############
#################################
"""
This script fit spline and regression difference to the data. This original fit will be used to create 
the generative noise model.
You can tune the spline parameters from the config file.
"""

def init_fit(config):

    savefig_path=pth.Path(config.simu_directory)    
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("simulation_kn"+str(knt))
        mkdir(savefig_path)        
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}
                
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                    kwargs_ml["mlfp"] = mlfp
                else:
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                saveml_path = saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                savefigml_path =saveml_path/"init_splines_and_resid/"

                #mkdir
                saveml_path = str(saveml_path)+"/"
                mkdir(saveml_path)
                savefigml_path = str(savefigml_path)+"/"
                mkdir(savefigml_path)
                
                lcs = config.get_lcs()
                
                if config.magshift is None :
                    magsft = [-np.median(lc.getmags()) for lc in lcs]
                else :
                    magsft = config.magshift
                pycs3.gen.lc_func.applyshifts(lcs, config.timeshifts, magsft) #remove median and set the time shift to the initial guess
                config.attachml(lcs, kwargs_ml)  # add microlensing
                
                spline = config.spl1(lcs, kn=knt)
                
                rls = pycs3.gen.stat.subtract(lcs, spline)
                
                chi_red = get_chi_red(spline,kwargs_ml,nl_lc=len(lcs))

                pycs3.gen.lc_func.display(lcs, [spline], showlegend=True, showdelays=True,filename=savefigml_path + "spline_fit.png")
                pycs3.gen.stat.plotresiduals([rls], filename=savefigml_path + "residual_fit.png")

                # and write data, again
                #if not os.path.isdir(config.lens_directory + config.combkw[i, j]):
                #    os.mkdir(config.lens_directory + config.combkw[i, j])

                pycs3.gen.util.writepickle((lcs, spline), saveml_path+'/initopt.pkl')
                
                with open(saveml_path+'chi_red.txt', 'w') as f:
                    f.write('chi_red: '+str(np.round(chi_red,4)))


# In[ ]:


################################
################################
################################
##### actually my analysis #####
################################
################################
################################

# > actually best to keep it separated

# Standard_analysis.ipynb


# In[5]:


#################################
########### Script 3a ###########
#################################

def run_DIC(lcs, spline, fit_vector, kn, mlparams, optim_directory, config, stream, tolerance=0.75):
    pycs3.sim.draw.saveresiduals(lcs, spline)
    print("I'll try to recover these parameters :", fit_vector)
    dic_opt = pycs3.pipe.optimiser.DicOptimiser(lcs, fit_vector, spline, config.attachml, mlparams, knotstep=kn,
                                    savedirectory=optim_directory,
                                    recompute_spline=True, max_core=config.max_core,
                                    n_curve_stat=config.n_curve_stat,
                                    shotnoise=config.shotnoise_type, tweakml_type=config.tweakml_type,
                                    tweakml_name=config.tweakml_name, display=config.display, verbose=False,
                                    correction_PS_residuals=True, max_iter=config.max_iter, tolerance=tolerance,
                                    theta_init=None)

    chain = dic_opt.optimise()
    dic_opt.analyse_plot_results()
    chi2, B_best = dic_opt.get_best_param()
    A = dic_opt.A_correction
    dic_opt.reset_report()
    dic_opt.report()

    if dic_opt.success:
        print("I succeeded finding a parameter falling in the %2.2f sigma from the original lightcurve." % tolerance)

    else:
        print("I didn't find a parameter that falls in the %2.2f sigma from the original lightcurve." % tolerance)
        print("I then choose the best one... but be carefull ! ")

    for k in range(len(lcs)):
        def tweakml_PS_NUMBER(lcs, spline):
            return twk.tweakml_PS(lcs, spline, B_PARAM, f_min=1 / 300.0, psplot=False, verbose=False,
                                  interpolation='linear', A_correction=A_PARAM)

        ut.write_func_append(tweakml_PS_NUMBER, stream,
                             B_PARAM=str(B_best[k][0]), NUMBER=str(k + 1), A_PARAM=str(A[k]))



# In[ ]:


def set_param_intrinsic(config):
    
    twk_dir = config.lens_directory+"figure/"
    mkdir(twk_dir)
    tweakml_plot_dir = twk_dir+'tweakml_plots/'
    mkdir(tweakml_plot_dir)    
    optim_directory = tweakml_plot_dir + 'twk_optim_%s_%s/' % (config.optimiser, config.tweakml_name)
    mkdir(optim_directory)

    savefig_path=pth.Path(config.simu_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("simulation_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                mllist_name,mlfp = ml_config 
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                    kwargs_ml["mlfp"] = mlfp
                else:
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                    
                combdir = saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                combdir = str(combdir)
                mkdir(combdir)
                
                if not os.path.isdir(combdir):
                    #os.mkdir(comb)
                    raise RuntimeError("This directory should already exist: ",combdir)
                    
                f = open(combdir + '/tweakml_' + config.tweakml_name + '.py', 'w+')
            
                f.write('import pycs3 \n')
                f.write('from pycs3.sim import twk as twk \n')    
            
                lcs, spline = pycs3.gen.util.readpickle(combdir+'/initopt.pkl' )
            
                pycs3.gen.mrg.colourise(lcs)
                pycs3.gen.lc_func.display(lcs,[spline],nicefont=True,showdelays=True,
                        filename=combdir+"/spline_fit.png")
                fit_vector = pycs3.pipe.optimiser.get_fit_vector(lcs, spline)  # we get the target parameter now

                # We need spline microlensing for tweaking the curve, if it is not the case we change it here to a flat spline that can be tweaked.
                # the resulting mock light curve will have no ML anyway, we will attach it the ML defined in your config file before optimisation.
                for k, l in enumerate(lcs):
                    if l.ml == None:
                        print('I dont have ml, I have to introduce minimal extrinsic variation to generate the mocks. Otherwise I have nothing to modulate.')
                        pycs3.gen.splml.addtolc(l, n=2)
                    elif l.ml.mltype == 'poly':
                        print('I have polyml and it can not be tweaked. I will replace it with a flat spline just for the mock light curve generation.')
                        l.rmml()

                # we replace the spline optimised with poly ml by one without ml 
                spline = pycs3.spl.topopt.opt_fine(lcs, nit=5, knotstep=knt,verbose=False, bokeps=knt / 3.0,
                                                       stabext=100) 
                for l in lcs:
                    pycs3.gen.splml.addtolc(l, n=2)
                pycs3.gen.util.writepickle((lcs, spline),combdir + '/initopt_generative_polyml.pkl')

                # Starting to write tweakml function depending on tweak_ml_type :
                run_DIC(lcs, spline, fit_vector, knt, kwargs_ml, optim_directory, config, f)
                list_string = 'tweakml_list = ['
                for k in range(len(lcs)):
                    list_string += 'tweakml_PS_' + str(k + 1) + ','
                list_string += ']'
                f.write('\n')
                f.write(list_string)
                f.close()
                
                ###
                ### These should only be relative to the reports which I don't consider
                ### for now ignore them
                ###
                """
                # rename the file :
                files = [file for file in os.listdir(optim_directory)
                         if os.path.isfile(os.path.join(optim_directory, file)) and (string_ML not in file)]
    
                for file in files:
                    prefix, extension = file.split('.')
                    os.rename(os.path.join(optim_directory, file),
                    os.path.join(optim_directory, prefix + "_kn%i_ml%s_%s%i." % (kn,ml[0], string_ML, ml[1]) + extension))
                """


# In[ ]:


#################################
########### Script 3b ###########
#################################

"""
This scrip will create copy of the data and mock light curves, according to your generative noise model.
I am using multithreading to do that.
"""
     
def draw_mock_para(config, combdir):
    current_dir = os.getcwd()
    os.chdir(combdir)
    lcs, spline = pycs3.gen.util.readpickle('initopt.pkl')

    pycs3.sim.draw.saveresiduals(lcs, spline)

    # add splml so that mytweakml will be applied by multidraw
    for l in lcs:
        if l.ml is None:
            pycs3.gen.splml.addtolc(l, n=2)
        elif l.ml.mltype == 'poly':
            polyml = True
            l.rmml()
            
    lcs, spline = pycs3.gen.util.readpickle('initopt_generative_polyml.pkl')
    pycs3.sim.draw.saveresiduals(lcs, spline)

    # import the module with the parameter of the noise :
    #print('I will use the parameter from : %s' % ('tweakml_' + config.tweakml_name + '.py'))
    exec(compile(open('tweakml_' + config.tweakml_name + '.py', "rb").read(),
                 'tweakml_' + config.tweakml_name + '.py', 'exec'), globals())

    multidraw(lcs, spline, onlycopy=False, n=config.nsim, npkl=config.nsimpkls,
                             simset=config.simset_mock, tweakml=tweakml_list,
                             shotnoise=config.shotnoise_type, trace=False,
                             truetsr=config.truetsr,trueMsr=config.trueMsr,
                             shotnoisefrac=1.0, scaletweakresi=False)
    os.chdir(current_dir)


def draw_mock_para_aux(arguments):
    return draw_mock_para(*arguments)

def my_draw_mocks(config):   
    if config.max_core is None:
        processes = cpu_count()
    else:
        processes = config.max_core

    p = Pool(processes=processes)
    print("Running on %i cores. " % processes)
    
    job_args = []
    
    savefig_path=pth.Path(config.simu_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("simulation_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                mllist_name,mlfp = ml_config 
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                else:
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                combdir = saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                combdir = str(combdir)

                simset = config.simset_mock
                file = glob.glob(os.path.join(combdir, "sims_" + simset + '/*.pkl'))
                if len(file) != 0:
                    for f in file:
                        os.remove(f)
                job_args.append((config,combdir))
    if processes > 1:
        p.map(draw_mock_para_aux, job_args)
    else:
        for args in job_args:
            draw_mock_para(*args)
    print("Done draw mocks")
    


# In[ ]:


#################################
########### Script 3c ###########
#################################


# I try to readapt the script from pycs 3c_optimise_copy_mocks.py so that I can use it even if the analysis is mine


"""
Optimise the copy and mock data. WARNING : this may take loooooooong. You probably want to launch that on several cores.
I'm not re-running on already optimized lcs ! It should be safe to launch this script many times, it will run on different batch of lightcurves.
"""


def exec_worker_mocks_aux(args):
    return exec_worker_mocks(*args)


"""def exec_worker_mocks(i, simset_mock, lcs, simoptfct, kwargs_optim, optset, tsrand, destpath,mltype,config):
    print("worker %i starting..." % i)
    time.sleep(i)
    sucess_dic = multirun(simset_mock, lcs, simoptfct, kwargs_optim=kwargs_optim,
                        optset=optset, tsrand=tsrand, keepopt=True, destpath=destpath,mltype=mltype,config=config)
    
    return sucess_dic
"""
def exec_worker_mocks(i, config, lcs,  kwargs_optim, optset, destpath):
    #j, config.simset_mock, lcs, config.simoptfct, kwargs, opts, config.tsrand, combdir
    print("worker %i starting..." % i)
    time.sleep(i)
    
    sucess_dic = multirun(config,lcs, kwargs_optim=kwargs_optim,
                        optset=optset,  keepopt=True, destpath=destpath)
    
    return sucess_dic

def write_report_optimisation(f, success_dic):
    if success_dic == None:
        f.write('This set was already optimised.\n')
    else:
        for i, dic in enumerate(success_dic):
            f.write('------------- \n')
            if dic == None:
                continue
            if dic['success']:
                f.write('None of the optimisations have failed for pickle %i. \n' % i)
            else:
                f.write('The optimisation of the following curves have failed in pickle %i : \n' % i)
                for id in dic['failed_id']:
                    f.write("   Curve %i :" % id + str(dic['error_list'][0]) + ' \n')
                f.write('\n')

def clean_previous_simopt(config, lcs,  kwargs, opts, combdir):
    # delete all previous optimisation
    destdir = os.path.join(combdir, "sims_%s_opt_%s" % (config.simset_mock, opts))
    try:
        listdir = os.listdir(destdir)
    except FileNotFoundError:
        return 0       
    for fi in listdir:
        os.remove(destdir+"/"+fi)
        
def optimise_sim(config):
    main_path = os.getcwd()
    base_lcs = config.get_lcs() 
    f = open(os.path.join(config.report_directory, 'report_optimisation_%s.txt' % (config.simoptfctkw)), 'w')

    savefig_path=pth.Path(config.simu_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("simulation_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                lcs = copy.deepcopy(base_lcs)
                mllist_name,mlfp = ml_config 
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                    kwargs_ml["mlfp"] = mlfp
                else:
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                    
                combdir = saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                combdir = str(combdir)
                if config.magshift is None :
                    magsft = [-np.median(lc.getmags()) for lc in lcs]
                else :
                    magsft = config.magshift

                timeshifts = config.timeshifts
                pycs3.gen.lc_func.applyshifts(lcs, timeshifts,magsft)  # be careful, this remove ml as well.
            
                # We also give them a microlensing model (here, similar to Courbin 2011)
                config.attachml(lcs, kwargs_ml)  # this is because they were saved as raw lcs, wihtout lcs.

                if config.max_core == None:
                    nworkers = cpu_count() -2
                else:
                    nworkers = config.max_core
                for c, opts in enumerate(config.optset):
                    kwargs = {'kn': knt, 'name': config.simoptfctkw}
                    kwargs = {**kwargs,**kwargs_ml} #combine the 2 dict
                    print("I will run the optimiser on the simulated lcs with the parameters :", kwargs)
                    p = Pool(nworkers)
                    #optfunc= lambda x: config.simoptfct(x)
                    #def optfunc(lcs,**kwargs):
                    #    return config.simoptfct(lcs,**kwargs)
                    #job_args = [(j, config.simset_mock, lcs, config.simoptfct, kwargs, opts, config.tsrand, combdir) for j in range(nworkers)]
                    job_args = [(j, config, lcs,  kwargs, opts, combdir) for j in range(nworkers)]
                    clean_previous_simopt(*job_args[0][1:])
                    success_list_simu = p.map(exec_worker_mocks_aux, job_args)
                    p.close()
                    p.join()
                    f.write('SIMULATIONS, kn%i, %s%i, optimiseur %s : \n' % (knt, ml_config[0],ml_config[1], kwargs['name']))
                    write_report_optimisation(f, success_list_simu)
                    f.write('################### \n')
    print("OPTIMISATION DONE : report written in %s" % (os.path.join(config.report_directory, 'report_optimisation_%s.txt' % config.simoptfctkw)))
    f.close()


# In[ ]:


#################################
########### Script 4a ###########
#################################
# plotting the error distribution mostly


# In[ ]:


#### From compare_sys.ipynb ########
def get_res_Group_i(config,knt,mltype,ml_config):
    """
    Get the result Group
    """
    
    savefig_path=pth.Path(config.analysis_directory)
    saveknt_path = savefig_path/str("analysis_kn"+str(knt))
            
    mllist_name,mlfp = ml_config 
    if mltype=="polyml":
        mlfp_str="_mlfp_"+str(mlfp)
    else:
        if config.forcen:
            mlfp_str="_nmlspl_"+str(mlfp)
        else:
            mlfp_str="_knstml_"+str(mlfp)
    """
    kwargs_ml = {"mltype":mltype,"mllist_name":mllist_name}
    if mltype=="polyml":
        mlfp_str="_mlfp_"+str(mlfp)
        kwargs_ml["mlfp"] = mlfp
    else:
        if config.forcen:
            mlfp_str="_nmlspl_"+str(mlfp)
        else:
            mlfp_str="_knstml_"+str(mlfp)
        kwargs_ml["forcen"] = config.forcen
        kwargs_ml["nmlspl"] = mlfp
    """
    data_path = saveknt_path/str("ml"+mltype[:-2]+mllist_name+mlfp_str)   
    data_path = str(data_path)
    sim_path  = data_path.replace("Analysis","Simulation").replace("analysis","simulation")
    sim_path = sim_path+"/sims_" + config.simset_mock+"_opt_"+config.optset[0]

    #ERROR
    error = Error(sim_path)
    error_distr = error.get_distr()
    error.create_error()
    #RES
    res_Group = getresults(data_path,error=error,labels=config.delay_labels)
    return res_Group


# In[ ]:


def combine_models(config,sigmathresh):
    wddir = config.combined_directory 
    marginalisation_plot_dir = wddir+'/figure/marginalisation_plots/'
    mkdir(marginalisation_plot_dir)        

    indiv_marg_dir = marginalisation_plot_dir + config.name_marg_spline + '/'
    mkdir(indiv_marg_dir)        

    marginalisation_dir = wddir+ config.name_marg_spline + '/'
    mkdir(marginalisation_dir)        

    
    f = open(marginalisation_dir + 'report_%s_sigma%2.1f.txt' % (config.name_marg_spline, sigmathresh), 'w')

    colors = config.colors
    color_id = 0

    group_list = []

    opt = config.optset[0]
        
    # Redefine the keyword here because you don't necessary want to marginalise over everything
    # mlknotsteps_marg now correspond to ml_config, which therefore contains all combinations of ml_lc and degrees

    masked_A = config.maskA
    
    savefig_path=pth.Path(config.analysis_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                mllist_name,mlfp = ml_config 
                group_i = get_res_Group_i(config,knt,mltype_i,ml_config)
                group_i.color=colors[color_id]
                name = config.combkw[knt_i][mlt_i][mlc_i].replace("_"," ")
                group_i.name = name
                group_list.append(group_i)
                color_id += 1
                if color_id >= len(colors):
                    color_id = 0  # reset the color form the beginning
                f.write('Set %s, knotstep : %2.2f, deg : %s %s \n' % (name, knt, mllist_name,mlfp))
                f.write('Tweak ml name : %s \n' % config.tweakml_name_marg_spline[0])
                f.write('------------------------------------------------ \n')
    
    #combine results
    
    #combined = combine_series(group_list, sigmathresh=sigmathresh)
    combined,comb_list,combined_indexes = combine_series_methodB(group_list, sigmathresh=sigmathresh,return_combined_list=True)
    
    
    print("Final combination for marginalisation ", config.name_marg_spline)

    savefig = indiv_marg_dir + config.name_marg_spline + "_sigma_%2.2f_myplot.pdf" % sigmathresh
    
    #create plot
    delayplot(group_list,savefig,colors=colors,refgroup=combined,selected_groups_indexes=combined_indexes)
    
    print("Saved group_list as ",str(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl'),\
        "and combined result as ",str(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl') )
    pkl.dump(group_list,
             open(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl',
                  'wb'))
    pkl.dump(combined,
             open(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl',
                  'wb'))
    
    #####################
    # combine mag shift #
    #####################
    mag_marginalisation_plot_dir = wddir+'/figure_MAG/marginalisation_plots/'
    mkdir(mag_marginalisation_plot_dir)
    
    mag_indiv_marg_dir = mag_marginalisation_plot_dir + config.name_marg_spline + '/'
    mkdir(mag_indiv_marg_dir)        
    mag_marginalisation_dir = wddir+ config.name_marg_spline + '/mag/'
    mkdir(mag_marginalisation_dir)   
    
    series_to_combine = []
    for name_to_comb in comb_list :
        str_to_comb = name_to_comb.replace(" ","_").replace("spl1_","") 
        
        knt = int(str_to_comb.split("_")[0].replace("ks",""))
        mltype =  str_to_comb.split("_")[1].replace("ml","")
        if "spl" in mltype:
            mlfp_ = "nmlspl_"
        elif "poly" in mltype:
            mlfp_ = "mlfp_"
        mlfp  = mlfp_+ str(str_to_comb.split("_")[-1])
        
        analysis_path = pth.Path(config.analysis_directory)
        anlknt_path  = analysis_path/str("analysis_kn"+str(knt))
        
        data_path = str(anlknt_path)+"/ml"+mltype+"_"+mlfp
        sim_path  = data_path.replace("Analysis","Simulation").replace("analysis","simulation")
        sim_path  = sim_path+"/sims_" + config.simset_mock+"_opt_"+config.optset[0]

        #ERROR
        error = Error_mag(sim_path)
        error.create_error()
        #RES
        res_Group = getresults_mag(data_path,error=error,labels=["AB","AC","BC"],name=name_to_comb)
        series_to_combine.append(res_Group)
        
    #mag_combined = series_to_combine[0]
    #if len(series_to_combine)>1:
    #    for G in series_to_combine[1:]:  
    #        mag_combined = combine_groups(mag_combined,G)   
    mag_combined = combine_group_list_methodB(series_to_combine)
    #mag_combined = combine_series_methodB(series_to_combine,sigmathresh=float("inf"))


    print("Final combination for marginalisation ", config.name_marg_spline)

    savefig = indiv_marg_dir +"mag_" +config.name_marg_spline + "_sigma_%2.2f_myplot.pdf" % sigmathresh
    
    #create plot
    dmagplot(series_to_combine,savefig,colors=colors,refgroup=mag_combined)
    
    print("Saved group_list as mag_",str(marginalisation_dir +  config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl'),\
        "and combined result as mag_",str(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl') )
    pkl.dump(series_to_combine,
             open(marginalisation_dir  +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl',
                  'wb'))
    pkl.dump(mag_combined,
             open(marginalisation_dir  +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl',
                  'wb'))



# In[ ]:


def standard_pipeline(config,verbose=False):    
    #Creating missing directories
    if verbose:
        print("Setting up directories\n")
    setup_directories(config)
    #Assuming Standard_analysis alread run, if not to implement a way to run it at this point
    # Simulation
    if verbose:
        print("Initial fit\n")
    init_fit(config)
    if verbose:
        print("Setting intrinsic parameter\n")
    set_param_intrinsic(config)
    if verbose:
        print("Drawing mocks\n")
    my_draw_mocks(config)
    # Optimise the simulation
    if verbose:
        print("Optimising the simulation\n")
    optimise_sim(config)
    # Plot resulting error distribution and original analysis distribution 
    if verbose:
        print("Plotting the error distribution\n")
    plot_err(config)
    
    print("Done standard error analysis\n")
    
    # Combine the results obtained from different models and plot them
    if verbose:
        print("Combining models\n")
    combine_models(config,sigmathresh=config.sigmathresh)


# In[ ]:


if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Find the noise parameter to reproduce the data. Then prepare the copies of the light curves and draw some mock curves. \
                               Analyse them. Obtain the error distribution, and the corresponding systematic. Repeat the process with the time delay corrected for the systematic",
                               formatter_class=ap.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process (Euler, SMARTS, ... )"
    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname
    verbose  = args.verbose

    present_program(sys.argv[0])


    config = get_config(lensname=lensname,dataname=dataname,config_path="myconfig")
    dt_string = time.strftime("%d/%m/%Y %H:%M:%S")
    print("Running: ",sys.argv[0])
    print("Machine: ",os.uname()[1]) 
    print("Config file:", lensname+"_"+dataname)
    print("Started the :", dt_string)
    print("#########################")
    standard_pipeline(config,verbose=verbose)

    #### 
    inspect_err_dir = config.combined_directory +'/figure/marginalisation_plots/'
    plt_err(config,savefig_dir=inspect_err_dir)
    plt_err_tot(config,savefig_dir=inspect_err_dir)
    plt_intr_err(config,savefig_dir=inspect_err_dir)
    #### 
    ######
    # FR mod
    #from observed_FR import combine_models_mag
    #combine_models_mag(config)
    ######
    print("\n",str(sys.argv[0]),": Done!")

