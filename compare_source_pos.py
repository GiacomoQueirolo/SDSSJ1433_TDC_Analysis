#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import numpy as np
import matplotlib.pyplot as plt

from tools import *
from plotting_tools import base_colors
from source_pos import get_source_pos_MCMC

if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Compare the reconstructed source pos in the source plane")
    #parser.add_argument("--no_corner_plot", action="store_false", dest="corner_plot", default=True,
    #                help="DO NOT plot the corner plot")
    parser.add_argument("-v","--verbose", action="store_true", dest="verbose", default=False,
                    help="Verbose")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings=args.SETTING_FILES
    #corner_plot = args.corner_plot
    verbose = args.verbose
    f,ax =plt.subplots(1, 1, figsize=(10,10))
    sets_names = ""
    sets_nm_short = ""
    for i,sets in enumerate(settings):
        color = base_colors[i]
        set_nm = strip_setting_name(sets)
        sets_names+=set_nm
        sets_nm_short+=set_nm+"_"
        if i!=len(settings):
            sets_names+=" and "
        if verbose:
            print(set_nm)
        kw_source,_ = get_source_pos_MCMC(sets,output_mcmc=False)
        src = kw_source["source_ra"][0],kw_source["source_dec"][0]
        src_errbound = list(kw_source["source_ra"][1:]),list(kw_source["source_dec"][1:])
        ax.scatter(*src,marker="x",label=set_nm+" source",c=color)
        ax.plot(src_errbound[0],src_errbound[1],color=color,linestyle="dashed")
        ax.plot(*src,marker="x",linestyle="dashed",c=color)
        ax.text(src[0],src[1],"Source",c=color)
    plt.title("Source position for "+sets_names)
    bck_dir = "backup_results/Source_sup/"
    mkdir(bck_dir)
    svnm = bck_dir+sets_nm_short[:-1]+".pdf"
    plt.savefig(svnm)
    print("Saving "+svnm)
    success()

