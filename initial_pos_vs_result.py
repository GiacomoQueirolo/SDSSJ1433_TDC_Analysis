#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Plot corner plot of posterior for position vs initial given positionf

import copy
import os, sys
import argparse
import numpy as np
from corner import corner 
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from get_res import *
from tools import *
from conversion import conv_radec_to_xy


# In[ ]:


if __name__=="__main__":
    ################################
    present_program(sys.argv[0])
    ################################
    parser = argparse.ArgumentParser(description="Corner plots of posterior of position of images, lens and perturber vs initial given pos. ")
    parser.add_argument("-pix","-p","-px", "--pixel", dest="pixel", default=False, action="store_true",
                        help="Give results in pixels")

    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args  = parser.parse_args()
    pix   = args.pixel
    setting_name =  args.SETTING_FILES
    
    backup_path = "backup_results"
    if pix:
        prm_to_analyse_nms = ["x main lens [pix]",   "y main lens [pix]",   "x prt [pix]",         "y prt [pix]"]
    else:
        prm_to_analyse_nms = ["ra main lens",        "dec main lens",       "ra prt",              "dec prt"]
    prm_to_analyse_no_sub  = ["center_x_lens_light0","center_y_lens_light0","center_x_lens_light1","center_y_lens_light1"]
    prm_to_analyse_sub     = ["center_x_lens0",      "center_y_lens0",      "center_x_lens_light0","center_y_lens_light0"]
    
    for sets in setting_name:
        seti     = get_setting_module(sets).setting()
        save_dir = get_savefigpath(sets)
        SUB      = check_if_SUB(sets) 
        
        if SUB:
            prm_names = prm_to_analyse_sub
        else:
            prm_names = prm_to_analyse_no_sub
        smpl = []
        for prm_name in prm_names:
            smpl.append(get_mcmc_smpl_for_prm(sets,prm_name,backup_path=backup_path))
        
        if pix:
            smpl_to_cnv = copy.deepcopy(smpl)
            smpl_cnv    = []
            # main_lens
            smpl_x,smpl_y         = conv_radec_to_xy(sets,ra=smpl_to_cnv[prm_to_analyse_sub.index("center_x_lens0")],\
                                                     dec=smpl_to_cnv[prm_to_analyse_sub.index("center_y_lens0")]) 
            # prt
            smpl_x_prt,smpl_y_prt = conv_radec_to_xy(sets,ra=smpl_to_cnv[prm_to_analyse_sub.index("center_x_lens_light0")],\
                                                     dec=smpl_to_cnv[prm_to_analyse_sub.index("center_y_lens_light0")]) 
            smpl = [smpl_x,smpl_y,smpl_x_prt,smpl_y_prt]
            
        init_val = []
        if not pix:
            init_val.append(seti.center_x)  #ra main lens
            init_val.append(seti.center_y)  #dec main lens
            init_val.append(seti.x_pert) #ra prt
            init_val.append(seti.y_pert) #dec prt
        
        else:
            ML_xy  = conv_radec_to_xy(sets,ra=seti.center_x, dec=seti.center_y)
            PRT_xy = conv_radec_to_xy(sets,ra=seti.x_pert,   dec=seti.y_pert)
            init_val.append(*ML_xy[0])  #x main lens
            init_val.append(*ML_xy[1])  #y main lens
            init_val.append(*PRT_xy[0]) #x prt
            init_val.append(*PRT_xy[1]) #y prt

        
        fg, ax = plt.subplots(len(prm_names),len(prm_names),figsize=(15,15))
        corner(np.array(smpl).T,labels=prm_to_analyse_nms,show_titles=True,truths=init_val,truth_color="r",fig=fg)
        print(init_val)
        init_val_str = str([np.round(i,2) for i in init_val]).replace("[","").replace("]","")
        legend_elements = [Patch(label="Initial value",facecolor="r"),Patch(label=init_val_str,facecolor="w")]
        axdel=ax[0][-1]
        axdel.legend(handles=legend_elements)
        axdel.axis("off")
        plt.tight_layout()
        if pix:
            print("Produced ",str(save_dir),"/compare_pos_results.pdf")
            fg.savefig(str(save_dir)+"/compare_pos_results.pdf")
        else:
            print("Produced ",str(save_dir),"/compare_pos_results_radec.pdf")
            fg.savefig(str(save_dir)+"/compare_pos_results_radec.pdf")
        #TODO: position of the images
    success(sys.argv[0])

