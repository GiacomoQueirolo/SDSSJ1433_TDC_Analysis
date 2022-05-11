import copy
import os, sys
import argparse
import warnings
import numpy as np
import pathlib as pth
from corner import corner
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#my libs
from tools import *
from get_res import *
from get_param_title import newProfile
from plotting_tools import my_corner_general,base_colors

parser = argparse.ArgumentParser(description="Plot the superposed corner plot of the posterior distribution of the Fermat potential difference and the other lens parameter from the given filters")
parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                    help="Cut the first <c> steps of the mcmc to ignore them")
parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the plot")
parser.add_argument("-if","--ignore_fermat",dest="ignore_fermat", default=False,action="store_true",
                    help="Ignore superposition of corner plot of fermat potential differences")
parser.add_argument("-ilp","--ignore_lens_param",dest="ignore_lens_param", default=False,action="store_true",
                    help="Ignore superposition of corner plot of lens paramaters")

parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

args = parser.parse_args()
cut_mcmc      = args.cut_mcmc
dir_name      = args.dir_name
setting_names = args.SETTING_FILES
ign_Df        = args.ignore_fermat
ign_lnsprm    = args.ignore_lens_param
if ign_Df and ign_lnsprm:
    print("No job to do == job done")
    exit()
    

setting_name = []
setting_position = []
for sn in setting_names:
    setting_position.append(find_setting_path(sn))
    setting_name.append(get_setting_name(sn))

  
#############################
present_program(sys.argv[0])
#############################


backup_path="backup_results"

filters  = [strip_setting_name(st) for st in setting_name]
save_dir = create_dir_name(setting_names,save_dir="PDF_Sup_Corner",dir_name=dir_name,backup_path=backup_path)

for st_i in range(len(setting_position)):
    os.system("cp "+str(setting_position[st_i])+"/"+str(get_setting_name(setting_names[st_i]))+    " "+str(save_dir)+"/.")
col = base_colors
if not ign_Df:
    print("Producing Df superposed corner plot")

    # for fermat potentials
    fermat_mcmc= [get_mcmc_fermat(st,backup_path) for st in setting_name]
    
    from fermat_pot_analysis import labels_Df as param_names
    param_names[-1] = param_names[0].replace("AB","BC")

    param_names = [ p+" [\"^2]" for p in param_names]
    
    samples_Df = []
    fg, ax = plt.subplots(3,3,figsize=(10,10))
    legend_elements  = []
    for i,mcmc_iT in enumerate(fermat_mcmc):
        cut_mcmc_scaled = int(len(mcmc_iT)*cut_mcmc/1000)
        mcmc_i   = np.transpose(mcmc_iT[cut_mcmc_scaled:]) 
        mcmc_Dfi = mcmc_i[1:]-mcmc_i[0]
        #############################################
        #############################################    
        #print("WARNING: Given the low S/N of image D, I will discard here and instead consider Delta BC")    
        #param_names=[r"$\Delta\phi AB$", r"$\Delta\phi AC$",r"$\Delta\phi AD$"]
        mcmc_c    = np.array(copy.deepcopy(mcmc_Dfi))
        mcmc_BC   = mcmc_c[1] - mcmc_c[0]  # BC = C - B = (C-A)-(B-A) = AC - AB
        mcmc_c[2] = mcmc_BC
        mcmc_Dfi  = mcmc_c 
        samples_Df.append(mcmc_Dfi.tolist())
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(setting_names[i])))
        if i==0:
            warnings.warn("Given the low S/N of image D, I will discard here and instead consider Delta BC")    
            corner(mcmc_Dfi.T,labels = param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
        else:
            corner(mcmc_Dfi.T,color=col[i%len(col)],fig=fg,plot_datapoints=False,hist_kwargs= {"density":True})
    axdel=ax[0][2]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    plt.tight_layout()
    fg.savefig(str(save_dir)+"/Df_sup.pdf")
    
    print("Produced Df_sup.pdf in "+str(save_dir))
    """
    #pt = my_corner(samples_Df,samples_names=setting_name,labels=param_names,udm="\"^2")
    pt = my_corner_general(samples_Df,samples_names=setting_name,labels=param_names)

    pt.savefig(str(save_dir)+"/superposed_Df.pdf")
    pt.close()
    print("Produced superposed_Df.pdf in "+str(save_dir))
    """
if not ign_lnsprm:
    print("Producing lens params superposed corner plot")

    smpl_mcmc = [get_mcmc_smpl(st,backup_path) for st in setting_name]
    prm_mcmc  = [get_mcmc_prm(st,backup_path) for st in setting_name]
    all_CP    = all([check_if_CP(st) for st in setting_name])

    param_to_compare= ['theta_E_lens0', 'e1_lens0', 'e2_lens0','theta_E_lens1', 'gamma_ext_lens2','psi_ext_lens2']
    if all_CP:
        param_to_compare = [param_to_compare[0],"gamma_lens0",*param_to_compare[1:] ]

    pair=[]
    for par in param_to_compare:
        pair_i=[]
        for i in range(len(prm_mcmc)):
            for j in range(len(prm_mcmc[i])):
                if par == prm_mcmc[i][j]:
                    pair_i.append(j)
                    break
        if len(pair_i)<len(prm_mcmc):
            warnings.warn("\nThis parameter is not present in all the settings (maybe fixed?): "+par+"\nIgnored for all settings\n")
            param_to_compare.remove(par)
            continue    

        pair.append(pair_i)


    samples_comp = [[] for i in setting_name]
    param_comp   = []
    prm_titles   = []

    for i in range(len(pair)):
        prm_titles.append(newProfile(param_to_compare[i])[0])
        for j in range(len(smpl_mcmc)):
            smpl_j = np.transpose(smpl_mcmc[j])[pair[i][j]]  # shape: steps
            if param_to_compare[i]=="psi_ext":
                samples_comp[j].append(smpl_j*180/np.pi)
            else:
                samples_comp[j].append(smpl_j)               # shape: setting, prm, steps
    fg, ax = plt.subplots(len(param_to_compare),len(param_to_compare),figsize=(15,15))
    legend_elements = []
    for i in range(len(samples_comp)):
        corner_data = np.array(samples_comp[i]).T 
        param_names = [newProfile(prm_to_cmp)[0] for prm_to_cmp  in param_to_compare]
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(setting_names[i])))
        if i==0:
            corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
        else:
            corner(corner_data,color=col[i%len(col)],fig=fg,plot_datapoints=False,hist_kwargs= {"density":True})
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    plt.tight_layout()
    fg.savefig(str(save_dir)+"/lnsprm_sup.pdf")
    
    print("Produced lnsprm_sup.pdf in "+str(save_dir))
    """
    pt = my_corner_general(samples_comp,samples_names=setting_name,labels=prm_titles)
    pt.savefig(str(save_dir)+"/superposed_lensprm.pdf")
    pt.close()

    print("Produced superposed_lensprm.pdf in "+str(save_dir))
    """
success(sys.argv[0])


