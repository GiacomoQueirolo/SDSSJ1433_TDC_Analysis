#!/usr/bin/env python
# coding: utf-8


# Redone for ABCD (also D)
# remodernisation done as well 

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
from stnd_red_chi import get_chi_red
from plot_distrib import plot_result_single
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
    
    orig_shift_time = config.timeshifts
    orig_shift_mag  = set_orig_magshift(lcs)
    
    savefig_path=pth.Path(config.analysis_directory)
    mkdir(savefig_path)
    
    
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        if verbose:
            print("analysing intrinsic knot ",knt)
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        mkdir(savefig_path)        
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            if verbose:
                print("Considering ml type ",mltype_i)

            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                if verbose:
                    print("Considering ml config ", ml_config)
                mllist_name,mlfp = ml_config 
                # note: mlfp is for polyml the n* of Free Parameters ofthe polynomial, or the 
                kwargs_mc = {"knst_list":config.knotstep, "nit":config.nit,
                 "dt_range":config.tsrand,"mc_res":config.mc_res}
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}

                if mltype_i=="polyml": 
                    kwargs_ml["mlfp"] = mlfp
                elif mltype_i=="splml": 
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                kwargs_mc.update(kwargs_ml)
                saveml_path      = saveknt_path/config.get_savemlpath(mltype_i,ml_config)# saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                savesplines_path = saveml_path/"splines/"

                #mkdir
                mkdir(saveml_path)
                saveml_path = str(saveml_path)+"/"
                mkdir(savesplines_path)
                savesplines_path = str(savesplines_path)+"/"
                print("Saving results in: "+saveml_path)

                
                with open(saveml_path+"/kwargs_mc.data","wb") as f:
                    pickle.dump(kwargs_mc,f) 
                data_vct = get_data_vect_single_analysis(config,savesplines_path,orig_shift_mag,knt,kwargs_ml)
                begin = time.time()
                
                if __name__ == '__main__':
                    with multiprocessing.Pool(n_processes) as p:
                        timedelays,dmags ,in_time_shift, in_mag_shift, chi_red, splines, residuals =\
                        zip(*p.map(partial(single_analysis,**data_vct),range(config.mc_res)),range(config.mc_res))
                end = time.time()
                #I don't know why but this might be necessary
                if timedelays[-1]==0:
                    if len(timedelays)==config.mc_res+1:
                        dmags  = dmags[:-1]
                        timedelays = timedelays[:-1]
                        in_mag_shift = in_mag_shift[:-1]
                        in_time_shift = in_time_shift[:-1]
                        splines = splines[:-1]
                        chi_red = chi_red[:-1]
                        residuals = residuals[:-1]
                if verbose: 
                    print("Done with multiprocessing in t=",(end-begin)/3600.," h\n")
    
                with open(saveml_path+"/td.data", 'wb') as f:
                    pickle.dump(timedelays, f)
                with open(saveml_path+"/dmags.data", 'wb') as f:
                    pickle.dump(dmags, f)
                with open(saveml_path+"/in_mag_shift.data", 'wb') as f:
                    pickle.dump(in_mag_shift, f)
                with open(saveml_path+"/in_time_shift.data", 'wb') as f:
                    pickle.dump(in_time_shift, f)
                with open(saveml_path+"/resid.data", 'wb') as f:
                    pickle.dump(residuals, f)
                with open(saveml_path+"/chi_red.data", 'wb') as f:
                    pickle.dump(chi_red, f)
                with open(saveml_path+"/splines.data", 'wb') as f:
                    pickle.dump(splines, f)
                # MOD_ROB
                success = check_analysis(timedelays,max_std=5) #max_std should be given by the config file
                with open(saveml_path+"/success.data","wb") as f:
                    pickle.dump(success,f)
                ########## Plot the resulting time delay distribution #############

                name = str(config.combkw[a][mlt_i][mlc_i]).replace("spl1_","")
                
                distr = np.transpose(timedelays)
                distr_dmag = np.transpose(dmags)
                plot_result_single(distr,name,saveml_path,labels=config.delay_labels)
                plot_result_single(distr_dmag,name,saveml_path,labels=config.delay_labels,mag=True)


#Auxiliary functions:
def set_orig_magshift(lcs):
    mag = []
    for i in range(len(lcs)):
        mg = lcs[i].getmags()
        mg_er = lcs[i].getmagerrs()
        mg_av = 0
        mg_err_sq=0
        for j in range( len(mg)):
            mg_av += mg[j]/(mg_er[j]**2)
            mg_err_sq+=1/(mg_er[j]**2)
        mag.append(mg_av/mg_err_sq) 
    true_shift_mag =-( np.array(mag) - mag[0])
    return true_shift_mag.tolist()

def mask_A_peak(lcs):
    lc=lcs[0]
    #MOD_MODERATE_MASK
    for i in [71,73,74,75, 83,84,85,86,87]:
        lc.magerrs[i]*=30
    return lcs
 

def get_data_vect_single_analysis(config,savesplines_path,orig_shift_mag,knt,kwargs_ml):
    orig_shift_time = config.timeshifts
    lcs = config.get_lcs() 
    if config.maskA:
        lcs=mask_A_peak(lcs)
    
    data_vect= {"spl":config.spl1,"mc_res":config.mc_res,"tsrand":config.tsrand,"lcs":lcs,\
            "get_mllist":config.get_mllist,\
            "attachml":config.attachml,\
            "savesplines_path":savesplines_path,\
            "orig_shift_mag":orig_shift_mag,\
            "orig_shift_time":orig_shift_time,"knt":knt,\
            "mlparams":kwargs_ml}

    return data_vect
#def single_analysis(config,mc_i,mc_res,tsrand,lcs,savesplines_path,orig_shift_mag,orig_shift_time,knt,mlparams,ret_lcs=False):
# change bc config is a module and cannot be pickled for multiprocess
def single_analysis(mc_i,spl,mc_res,tsrand,lcs,get_mllist,attachml,savesplines_path,orig_shift_mag,orig_shift_time,knt,mlparams,ret_lcs=False):
    if mc_i%(20/mc_res)==0:
        print("Iteration "+str(mc_i+1)+" of "+str(mc_res)+" ("+str(np.round((mc_i+1)*100/mc_res,3))+"%)")
    
    np.random.seed(mc_i)
    
    sig_mag =[np.median(lci.getmagerrs()) for lci in lcs]

    #this avoid that bad choices of parameters make the spline fitting fail
    while True:
        # initialise mag, dt and lightcurves:
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
        
        # consider the microlensing for lc B, C and D
        attachml(lcs,mlparams) 
        try: 
            spline = spl(lcs,kn=knt)
            break
        except:
            pass
         
    timedelays  = get_dts(lcs,only_indep=True,wD=True)
    # Auxiliary info # discontinued for now
    # so we can see the spline and ml when secondary peaks appears
    """
    if np.any([np.abs(list(timedelays.values())[k] + orig_shift_time[1:][k])>7 for k in range(len(orig_shift_time))]):
        file_lcs_bad =str(savesplines_path)+"/large_ddt_"+str(mc_i+1)+".data"
        with open(file_lcs_bad,"wb") as f:
            pickle.dump([*lcs,spline,all_dt,mlparams],f)
    """

    #initial time and mag shifts
    in_time_shift = shift_time
    in_mag_shift  = [lc.magshift for lc in lcs]
    
    # FR study - 16th March
    mllist = get_mllist(mlparams["mllist_name"] )
    mag_ml_shift = []
    for j in range(len(lcs)):
        if j in mllist:
            if mlparams["mltype"]=="polyml":
                mag_ml_shift.append(lcs[j].ml.getfreeparams()[-1])
            elif mlparams["mltype"]=="splml":
                # if spline, I consider the average of the spline as a mean shift
                lc = lcs[j]
                dt=lc.jds-lc.jds[0]
                points=lc.ml.spline.eval(dt) 
                mag_ml_shift.append(np.mean(points))
        else:
            mag_ml_shift.append(0)
    mag_shift = np.array(mag_ml_shift) + np.array(in_mag_shift) #to be ADDED to the mag of the lcs in order to shift them
    # mag_shift is how much the lc must be shifted to superpose with the intrinsic lc
    # hence it is the inverser (x-1) of the magnitude shift operated by the lens
    # which we want to measure-> we call it dmag st: dmag=-mag_shift
    dmag = -mag_shift
    # (Dmag will be the differences of dmag btw images: Dmag_i = dmag_i - dmag_A)    
    #raise # check and decide how to do the next step    
    Dmag = dmag[1:]-dmag[0] # AB,AC,AD
    Dmag = Dmag.tolist()
    
    #residuals and chired
    residual = pycs3.gen.stat.subtract(lcs,spline)
    chi_red  = get_chi_red(spline,mlparams,nl_lc=len(lcs))
    # we save an example of spline
    if mc_i==0: 
        pycs3.gen.lc_func.display(lcs,[spline],nicefont=True,showdelays=True,
                    filename=savesplines_path+"fit.png")
    
    if ret_lcs:
        return lcs,timedelays,Dmag,\
        in_time_shift,in_mag_shift, chi_red, spline, residual
    else:
        return timedelays,Dmag,in_time_shift,in_mag_shift,\
            chi_red, spline, residual

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
    dts_lett = list(kwdt.keys())
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
    ####
    #comand = "python observed_FR.py "+lensname+" "+dataname+" &> logs/log_FR_"+lensname+"_"+dataname+".log &" 
    #os.system(comand)
    ####

