#!/usr/bin/env python
# coding: utf-8

# In[2]:

import sys,os
import numpy as np
import argparse as ap
import pathlib as pth

from Utils.tools import *
from TD_analysis.stnd_handling_data import getresults, Error


def check_robustness(config,verbose,threshold=5):
    savefig_path=pth.Path(config.analysis_directory)

    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                # note: mlfp is for polyml the n* of Free Parameters ofthe polynomial, or the 
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                else:
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                        
                saveml_path = saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)

                
                timedelays = getresults(saveml_path).data
                std_timedelays = np.std(timedelays,axis=1)
                

                name_ml = " ".join(str(saveml_path).split("/")[3:]).replace("analysis_","").replace("ml","").replace("_"," ")
 
                
                sim_path  = str(saveml_path).replace("Analysis","Simulation").replace("analysis","simulation")
                sim_path  = sim_path+"/sims_" + config.simset_mock+"_opt_"+config.optset[0]

                error = Error(sim_path)
                error_distr = error.get_distr()
                error.create_error()
                if verbose:
                    print(name_ml)
                    print("intrinsic variation: ",np.round(std_timedelays,3))
                    print("total error:         ",np.round(error.tot,3))
                for intrinsic,general in zip(std_timedelays,error.tot):
                    ratio = general/intrinsic
                    if verbose:
                        print("Ratio: ~",int(np.round(ratio)))
                    if ratio<threshold:
                        print(name_ml," has intrinsic error too large")
                        raise
                if verbose:
                    print()
    print("All successful!")

if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Robustness check: Check scatter of analysed time delays and compare it with uncertainties",
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

    config =get_config(lensname,dataname)
    check_robustness(config,verbose)
    success(sys.argv[0])

