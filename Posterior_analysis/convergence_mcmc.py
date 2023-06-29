# following Etl23:
# we have a small convergence definition
# equ. 7
import numpy as np
from Utils.get_res import get_mcmc_logL
from Utils.tools import *

import sys
import argparse

def test_convergence(setting,mcmc_logL=None,crit=5,crit_strict=2):
    if mcmc_logL is None:
        mcmc_logL     = get_mcmc_logL(setting)
    #cut_logL = mcmc_logL[int(len(logL)/2.):]
    cut_index=int((2000/100000)*len(mcmc_logL))
    med1 = np.median(mcmc_logL[:cut_index])
    med2 = np.median(mcmc_logL[-cut_index:])
    DlogL = abs(med2-med1)
    print(f"DlogL = {DlogL}")
    if DlogL<=crit:
        print("CONVERGED!")
        if DlogL<=crit_strict:
            print("Passed even the strict convergence criterium!")
    else:
        print("not converged :(")

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Check convergence criterium for MCMC")
    parser.add_argument('SETTING_FILES',nargs="+", help="Setting files to consider")
    present_program(sys.argv[0])

    args     = parser.parse_args()
    setting  = get_setting_module(args.SETTING_FILES,1)
    for sett in setting:
        print_setting(sett)
        test_convergence(setting=sett)
    success(sys.argv[0])    
