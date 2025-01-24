#!/usr/bin/env python
# coding: utf-8


# same as Standard_analysis_ABCD, but done for the gaussian kernel regression difference

import time
import pickle
import os,sys
import importlib
import numpy as np
import pathlib as pth
import argparse as ap
from copy import copy
import multiprocessing  
from functools import partial
import matplotlib.pyplot as plt
from numpy.random import normal as norm
from numpy.random import uniform as unif
# mod for 2nd chance
from astropy.stats import sigma_clip

import pycs3.gen.mrg
import pycs3.gen.stat
import pycs3.gen.splml
import pycs3.spl.topopt
import pycs3.gen.lc_func

########################
from Utils.tools import *
from stnd_red_chi import get_chi_regdiff
from plot_distrib import plot_result_single,plot_regdiff_fit
########################




# to have standardised plots
plt.rcParams['font.size']= 16
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18


n_processes= multiprocessing.cpu_count()


#MOD_ROB
def check_analysis(td_data,max_std=1.1,second_chance=False):
    # max_std = maximum acceptable standard deviation [d]
    distr = np.transpose(td_data)
    for i in range(len(distr)): #ABCD          
        ddt = distr[i]
        std_i = np.std(ddt)
        if std_i>max_std:
            # image D have a low prob second peak that 
            # draw all result to be discarded
            # in order to still consider that, we consider a 2nd chance by 
            # considering a sigma-clipped analysis
            if second_chance:
                clipped_ddt = sigma_clip(ddt, maxiters=5,copy=False,masked=False)
                std_i = np.std(clipped_ddt)
                if std_i>max_std:
                    return {"success":False,"std_i":std_i,"lc":i,"max_std":max_std,"second_chance":second_chance}
                else: 
                    continue
            return {"success":False,"std_i":std_i,"lc":i,"max_std":max_std,"second_chance":second_chance}
    return  {"success":True,"max_std":max_std,"second_chance":second_chance}
                    
def standard_analysis(config,verbose=False):
    config = get_config(config)
    
    lcs = config.get_lcs() 
    
    if config.maskA:
        lcs=mask_A_peak(lcs)
        
    if verbose:
       # Useful info
        lc_i = lcs[0]
        print ( "Lenght of observation campaign: ",(max(lc_i.getjds())-min(lc_i.getjds()))," days")
        print ( "N* of nights of observation campaign: ",len(lc_i.getjds()) ," days")
    
    orig_shift_mag  = set_orig_magshift(lcs)
    savefig_path    = pth.Path(config.regdiff_directory)
    mkdir(savefig_path)
    for rgd_set in config.regdiff_sets:
        saveset_path = savefig_path/rgd_set.set_name
        mkdir(saveset_path)
        savefit_path = saveset_path/"regdiff_fit/"
        mkdir(savefit_path)
        saveset_path = str(saveset_path)
        print("Saving results in: ",saveset_path)
        
        kwargs_mc = rgd_set.get_kw()
        kwargs_mc["orig_shift_mag"]  = orig_shift_mag 
        kwargs_mc["orig_shift_time"] = config.timeshifts 
        kwargs_mc["mc_res"]          = config.mc_res
        kwargs_mc["tsrand"]          = config.tsrand
        kwargs_mc["regdiff"]         = config.regdiff1
        kwargs_mc["savefit_path"]    = savefit_path

        name = kwargs_mc["set_name"]
        del kwargs_mc["set_name"]
        with open(saveset_path+"/kwargs_mc.data","wb") as f:
            pickle.dump(kwargs_mc,f) 

        begin = time.time()
        if __name__ == '__main__':
            with multiprocessing.Pool(n_processes) as p:
                #timedelays,dmags ,in_time_shift, in_mag_shift, chi_red, rslcs, residuals =\
                timedelays,dmags ,in_time_shift, in_mag_shift, chi, rslcs =\
                    zip(*p.map(partial(single_analysis_regdiff,**kwargs_mc),range(config.mc_res)),range(config.mc_res))
        end = time.time()
        if verbose: 
            print("Done with multiprocessing in t=",(end-begin)/3600.," h\n")
        #I don't know why but this might be necessary
        if timedelays[-1]==0:
            if len(timedelays)==config.mc_res+1:
                dmags         = dmags[:-1]
                timedelays    = timedelays[:-1]
                in_mag_shift  = in_mag_shift[:-1]
                in_time_shift = in_time_shift[:-1]
                chi           = chi[:-1]
                rslcs         = rslcs[:-1]

        with open(saveset_path+"/td.data", 'wb') as f:
            pickle.dump(timedelays, f)
        with open(saveset_path+"/dmags.data", 'wb') as f:
            pickle.dump(dmags, f)
        with open(saveset_path+"/in_mag_shift.data", 'wb') as f:
            pickle.dump(in_mag_shift, f)
        with open(saveset_path+"/in_time_shift.data", 'wb') as f:
            pickle.dump(in_time_shift, f)
        #with open(saveml_path+"/resid.data", 'wb') as f:
        #    pickle.dump(residuals, f)
        with open(saveset_path+"/chi.data", 'wb') as f:
            pickle.dump(chi, f)
        with open(saveset_path+"/rslcs.data", 'wb') as f:
            pickle.dump(rslcs, f)
        # MOD_ROB
        success = check_analysis(timedelays,max_std=5) #max_std should be given by the config file
        with open(saveset_path+"/success.data","wb") as f:
            pickle.dump(success,f)
        ########## Plot the resulting time delay distribution #############

        
        distr      = np.transpose(timedelays)
        distr_dmag = np.transpose(dmags)
        plot_result_single(distr,name,saveset_path,labels=config.delay_labels,rgd_set=rgd_set)
        plot_result_single(distr_dmag,name,saveset_path,labels=config.delay_labels,mag=True,rgd_set=rgd_set)


#Auxiliary functions:
def set_orig_magshift(lcs):
    mag = []
    for i in range(len(lcs)):
        mg        = lcs[i].getmags()
        mg_er     = lcs[i].getmagerrs()
        mg_av     = 0
        mg_err_sq = 0
        for j in range( len(mg)):
            mg_av     += mg[j]/(mg_er[j]**2)
            mg_err_sq += 1/(mg_er[j]**2)
        mag.append(mg_av/mg_err_sq) 
    true_shift_mag =-( np.array(mag) - mag[0])
    return true_shift_mag.tolist()

def mask_A_peak(lcs):
    lc=lcs[0]
    #MOD_MODERATE_MASK
    for i in [71,73,74,75, 83,84,85,86,87]:
        lc.magerrs[i]*=30
    return lcs
 
def single_analysis_regdiff(mc_i,mc_res,lcs,regdiff,tsrand,orig_shift_mag,orig_shift_time,savefit_path,ret_lcs=False,**kwargs_reg):
    if mc_i%(20/mc_res)==0:
        print("Iteration "+str(mc_i+1)+" of "+str(mc_res)+" ("+str(np.round((mc_i+1)*100/mc_res,3))+"%)")
    np.random.seed(mc_i)
    sig_mag =[np.median(lci.getmagerrs()) for lci in lcs]

    shift_mag  = copy(orig_shift_mag)    
    shift_time = copy(orig_shift_time)  
    for l in lcs:
        l.resetshifts() 
    # random variation of initial guesses
    for j in range(len(lcs)):
        dt_unif = unif(-tsrand,tsrand)
        shift_time[j] = shift_time[j] + dt_unif

        dmag_norm    = norm(0,sig_mag[j] )
        shift_mag[j] = shift_mag[j] + dmag_norm 

    # we had an initial guess of the time shift and mag shift
    pycs3.gen.lc_func.applyshifts(lcs, shift_time, shift_mag)

    # computation
    fits_rslcs, error_fct = regdiff(lcs, **kwargs_reg)
    timedelays = get_dts(lcs=lcs,only_indep=True,wD=True)
    #initial time and mag shifts
    in_time_shift = shift_time
    in_mag_shift  = [lc.magshift-lcs[0].magshift for lc in lcs[1:]] # i - A 
    # FR study, adapted for RegDiff
    rslcs_residuals = [pycs3.regdiff.rslc.subtract(fits_rslcs[i],fits_rslcs[0]) for i in range(1,len(fits_rslcs))] # i - A
    reg_mag_shift   = [rslcs_residuals[i].mags.mean() for i in range(len(rslcs_residuals))]

    # I HAVE TO DOUBLE CHECK THE MAG COMPUTATION, I AM NOT SURE OF IT
    dmag = np.array(reg_mag_shift) + np.array(in_mag_shift) #to be ADDED to the mag of the lcs in order to shift them
    # dmag HERE is the difference of mag btw the fits for lc I wrt lc A
    # hence it is the inverser (x-1) of the magnitude shift operated by the lens
    # which we want to measure-> we call it dmag st: dmag=-mag_shift
    #dmag = -mag_shift
    Dmag = dmag.tolist()
    
    #chi : not reduced (yet)
    chi  = np.mean([get_chi_regdiff(lc,rslc) for lc, rslc in zip(lcs,fits_rslcs)])
    # we save an example of spline
    if mc_i==0: 
        plot_regdiff_fit(lcs,fits_rslcs,rslcs_residuals,name=savefit_path+"/regdiff.pdf")

    if ret_lcs:
        results = lcs,timedelays,Dmag ,in_time_shift, in_mag_shift, chi, fits_rslcs 
    else:
        results = timedelays,Dmag ,in_time_shift, in_mag_shift, chi, fits_rslcs 
    return results


def get_kwdt(lcs):
    dts_string = pycs3.gen.lc_func.getnicetimedelays(lcs)
    kwdt       = {}
    for dts in dts_string.split("\n"):
        IJ  = dts[:2]
        dts = dts.split("=")[1]
        kwdt[IJ] = -float(dts)
    return kwdt

def get_dts(lcs,only_indep=False,wD=True):
    kwdt = get_kwdt(lcs)
    dts  = np.array(list(kwdt.values()))
    if not only_indep:
        return dts
    #dts_lett = list(kwdt.keys())
    # we consider A as the reference
    #index_A = np.array([ "A" in ksi for ksi in dts_lett])
    # better to just give it hardcoded
    if wD:
        dts = [kwdt["AB"],kwdt["AC"],kwdt["AD"]]
    else:
        dts = [kwdt["AB"],kwdt["AC"],kwdt["BC"]]
    return dts
    
if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="My analysis",
                               formatter_class=ap.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process"
    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args      = parser.parse_args()
    lensname  = args.lensname
    dataname  = args.dataname
    verbose   = args.verbose
    name_prog = sys.argv[0]
    present_program(name_prog)

    sys.path.append("myconfig/")
    config_file = "myconfig_" + lensname + "_" + dataname
    config = importlib.import_module(config_file)
    mkdir(config.lens_directory)
    
    dt_string = time.strftime("%d/%m/%Y %H:%M:%S")
    print("Machine: ",os.uname()[1]) 
    print("Config file:", config_file)
    print("Started the :", dt_string)
    print("#########################")
    standard_analysis(config,verbose=verbose)


