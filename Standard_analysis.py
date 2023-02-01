#!/usr/bin/env python
# coding: utf-8
# Mod from Notes 19th october: we consider a spline ml and see the difference with the only polynomial case 
#-> Expected: larger error distribution, smaller systematic
# In[1]:


import time
import os,sys
import importlib
import json,pickle
import numpy as np
import pathlib as pth
import argparse as ap
from copy import copy
import multiprocessing 
from corner import quantile
from datetime import datetime
from functools import partial
import matplotlib.pyplot as plt
from numpy.random import normal as norm
from numpy.random import uniform as unif

import pycs3.gen.mrg
import pycs3.gen.stat
import pycs3.gen.splml
import pycs3.spl.topopt
import pycs3.gen.lc_func

########################
from tools import *
from stnd_red_chi import get_chi_red
from plot_distrib import plot_result_single
########################


# In[ ]:


# to have standardised plots
plt.rcParams['font.size']= 16
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18


# In[ ]:


n_processes= multiprocessing.cpu_count()


# In[6]:


#MOD_ROB
def check_analysis(td_data,max_std=1.1):
    # max_std = maximum acceptable standard deviation [d]
    distr = np.transpose(td_data)
    for i in range(len(distr)): #ABCD          
        ddt = distr[i]
        std_i = np.std(ddt)
        if std_i>max_std:
            return {"success":False,"std_i":std_i,"lc":i}
    return  {"success":True}
                    
def standard_analysis(config,verbose=False):
    
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
                    mlfp_str="_mlfp_"+str(mlfp)
                    kwargs_ml["mlfp"] = mlfp
                elif mltype_i=="splml":
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
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
                    
                data_vct= {"mc_res":config.mc_res,"tsrand":config.tsrand,"lcs":lcs,\
                           "savesplines_path":savesplines_path,\
                           "orig_shift_mag":orig_shift_mag,\
                           "orig_shift_time":orig_shift_time,"knt":knt,\
                           "mlparams":kwargs_ml}
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
                success = check_analysis(td_data)
                with open(saveml_path+"/success.data","wb") as f:
                    pickle.dump(success,f)
                ########## Plot the resulting time delay distribution #############
                name = str(config.combkw[a][mlt_i][mlc_i]).replace("spl1_","")
                distr = np.transpose(timedelays)
                distr_dmag = np.transpose(dmags)
                plot_result_single(distr,name,saveml_path)
                plot_result_single(distr_dmag,name,saveml_path,mag=True)


# In[ ]:


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

#@jit
def single_analysis(mc_i,mc_res,tsrand,lcs,savesplines_path,orig_shift_mag,orig_shift_time,knt,mlparams):
    if mc_i%(20/mc_res)==0:
        print("Iteration "+str(mc_i+1)+" of "+str(mc_res)+" ("+str(np.round((mc_i+1)*100/mc_res,3))+"%)")
    
    np.random.seed(mc_i)
    
    sig_mag =[np.median(lci.getmagerrs()) for lci in lcs]

    
    spl = config.spl1
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
        config.attachml(lcs,mlparams) 
        try: 
            spline = spl(lcs,kn=knt)
            break
        except:
            pass
        
    l = pycs3.gen.lc_func.getnicetimedelays(lcs) 
    all_dt=[]
    for j in (l.split(" ")):
        if '\n' in j or j==l.split(" ")[-1]:
            all_dt.append(-float(j.split("\n")[0]))  
    timedelay_ab,timedelay_ac,timedelay_bc = all_dt[:3]
    
    # Auxiliary info
    # so we can see the spline and ml when secondary peaks appears
    if np.any([np.abs(all_dt[k] +orig_shift_time[1:][k])>7 for k in range(len(all_dt[:2]))]):
        str_dt = str(np.round(all_dt,0)).replace(" ","").replace("[","").replace("]","").replace(".","_")[:-1]
        file_lcs_bad =str(savesplines_path)+"/large_ddt_"+str(mc_i+1)+".data"
        with open(file_lcs_bad,"wb") as f:
            pickle.dump([*lcs,spline,all_dt,mlparams],f)
            
    #initial time and mag shifts
    in_time_shift = shift_time
    in_mag_shift  = [lc.magshift for lc in lcs]
    
    # FR study - 16th March
    mllist = config.get_mllist(mlparams["mllist_name"] )
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
    # which we want to measure-> we call it dmag st: dmag=-magshift
    # (Dmag will be the differences of dmag btw images: Dmag_i = dmag_i - dmag_A)
    dmag = -mag_shift
    dmag = dmag.tolist()

    
    #residuals and chired
    residual = pycs3.gen.stat.subtract(lcs,spline)
    chi_red  = get_chi_red(spline,mlparams,nl_lc=len(lcs))
    # we save an example of spline
    if mc_i==0: 
        pycs3.gen.lc_func.display(lcs,[spline],nicefont=True,showdelays=True,
                    filename=savesplines_path+"fit.png")
    
    return [timedelay_ab,timedelay_ac,timedelay_bc],dmag,\
            in_time_shift,in_mag_shift, chi_red, spline, residual


# In[ ]:


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

    success(name_prog)


# In[ ]:




