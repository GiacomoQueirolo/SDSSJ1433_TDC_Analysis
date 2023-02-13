import argparse


import matplotlib.pyplot as plt
import json
from numpy import inf
import numpy as np
import sys

from Utils.tools import *
from Utils.get_res import *

def plot_PSO(setting_name,backup_path="./backup_results/",verbose=False,logL=True):
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    savefig_path  = get_savefigpath(setting_name,backup_path)
    distance,likelihood= load_whatever(savemcmc_path+"/mypsosampling.json")


    if all(l==-inf for l in likelihood):
        raise RuntimeError("All likelihood are -inf")
        
    l = []
    d = []
    for i,dis in zip(likelihood,distance):
        if i != -inf:
            l.append(i)
            d.append(dis)
    if logL : 
       l = np.log(-np.array(l))
    #min_lkl  = np.min(l)
    if logL:
        res      = np.min(l)
        logL_fac = +1
        rnd      = 3
    else:
        res      = np.max(l)
        logL_fac = -1
        rnd      = 0
    #diffmm   = res - min_lkl
    min_lkl  = res+0.01*abs(res)*logL_fac
    dist_res = d[np.where(l==res)[0][0]]
    max_d    = d[np.where(abs(l-min_lkl)==np.min(abs(l-min_lkl)))[0][0]]       
    #likelihood = [min_lkl if x==-inf else x for x in likelihood]
    #plt.ylim(res-.03*diffmm,res+.001*diffmm)
    plt.figure(figsize=(12,9)) 
    if logL:
        plt.ylim(res-0.0001*abs(res),min_lkl)
        plt.xlim(min(d)-0.01*max_d,max_d)
    else:
        plt.ylim(min_lkl,res+0.0001*abs(res))
        plt.xlim(min(d)-0.01*max_d,max_d)
    plt.scatter(d,l,alpha=.1,marker=".",label="Optimum="+str(np.round(res,rnd))+" found at d="+str(dist_res))
    plt.axhline(res,c="r")
    plt.legend()
    plt.xlabel("Normalised distance from PSO result []")
    if logL:
        plt.ylabel("Log(Likelihood) []")
        plt.title("LogL vs Distance to optimum for best PSO points")
    else:
        plt.ylabel("-Likelihood []")
        plt.title("-Likelihood vs Distance to optimum for best PSO points")
    plt.savefig(savefig_path+"/psos_plot.png")
    if verbose:
        print("PSO plot saved as "+savefig_path+"/psos_plot.png")
    

if __name__=="__main__":
    present_program(sys.argv[0])
    parser = argparse.ArgumentParser(description="Produce, if the data is present, the plot of the PSO fitnesses with respect to the distance of the particles from the final result")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    parser.add_argument('-l','--Like',default=True,dest="logL",action="store_false",help="Consider -likelihood instead of log(-like)")
    parser.add_argument('-v','--verbose',default=True,dest="verbose",action="store_false",help="Not verbose")
    args = parser.parse_args()
    logL = args.logL
    verbose = args.verbose
    settings=args.SETTING_FILES
    for sets in settings:
        if verbose:
            if logL:
                print("Considering logL ...")
            else:
                print("Considering -Likelihood ...")
        plot_PSO(sets,verbose=verbose,logL=logL)
    success(sys.argv[0])
