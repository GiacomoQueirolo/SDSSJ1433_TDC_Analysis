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
from Utils.tools import *
from Utils.get_res import *
from Utils.get_param_title import newProfile

from Data.Param import get_prm_list
from Data.conversion import qphi_from_e1e2

from Plots.plotting_tools import my_corner_general,base_colors
from Posterior_analysis.fermat_pot_analysis import labels_Df as  param_names    
#################################################################
fnt = 16
plt.rcParams['xtick.labelsize'] = fnt
plt.rcParams['ytick.labelsize'] = fnt 
plt.rcParams['font.size'] = fnt
plt.rc('axes', labelsize=fnt)    # fontsize of the x and y labels
plt.rc('font', size=fnt)          # controls default text sizes
plt.rc('axes', titlesize=fnt)     # fontsize of the axes title
plt.rc('axes', labelsize=fnt)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=fnt)    # fontsize of the tick labels
plt.rc('ytick', labelsize=fnt)    # fontsize of the tick labels
plt.rc('legend', fontsize=fnt)    # legend fontsize
#################################################################
    
def plot_sup_corner_Df(setting_list,fermat_mcmc=None,col=base_colors,savefig_dir=None,name_pdf="Df_sup",cut_mcmc=None,simple_legend=False,
        param_names=param_names,already_BC=False,backup_path="backup_results/"):
    print("Producing Df superposed corner plot")
    # for fermat potentials
    settings_name = get_setting_name(setting_list)
    if fermat_mcmc is None:
        fermat_mcmc = [get_mcmc_fermat(st,backup_path) for st in settings_name]
        already_BC  = False
    param_names[-1] = param_names[0].replace("AB","BC")
    param_names     = [ p.replace(" [\"^2]","")+" [\"^2]" for p in param_names]
    
    samples_Df = []
    fg, ax = plt.subplots(3,3,figsize=(10,10))
    legend_elements  = []
    for i,mcmc_iT in enumerate(fermat_mcmc):
        cut_mcmc_scaled = 0
        if cut_mcmc:
            cut_mcmc_scaled = int(len(mcmc_iT)*cut_mcmc/1000)
        mcmc_i   = np.transpose(mcmc_iT[cut_mcmc_scaled:]) 
        if not already_BC:
            mcmc_Dfi = mcmc_i[1:]-mcmc_i[0]
            #############################################
            #############################################    
            mcmc_c    = np.array(copy.deepcopy(mcmc_Dfi))
            mcmc_BC   = mcmc_c[1] - mcmc_c[0]  # BC = C - B = (C-A)-(B-A) = AC - AB
            mcmc_c[2] = mcmc_BC
            mcmc_Dfi  = mcmc_c 
        else:
            mcmc_Dfi  = mcmc_i
        samples_Df.append(mcmc_Dfi.tolist())
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(settings_name[i],filter=simple_legend)))
        if i==0:
            warnings.warn("Given the low S/N of image D, I will discard here and instead consider Delta BC")    
            corner(mcmc_Dfi.T,labels = param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
        else:
            corner(mcmc_Dfi.T,color=col[i%len(col)],fig=fg,plot_datapoints=False,hist_kwargs= {"density":True})
    axdel=ax[0][2]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    #plt.tight_layout()
    if savefig_dir:
        fg.savefig(str(savefig_dir)+"/"+name_pdf+".pdf")
        print("Produced "+name_pdf+".pdf in "+str(savefig_dir))
    else:
        return fg,ax

def plt_sup_corner_lnsprm(setting_list,smpl_mcmc=None,prm_mcmc=None,cut_mcmc=None,stnd_lnsprm=False,ign_lnsprm_qphi=False,
                    col=base_colors,savefig_dir=None,simple_legend=False,name_pdf="lnsprm_sup",backup_path="backup_results/"):
    print("Producing lens params superposed corner plot")
    setting_name = get_setting_name(setting_list) 
    samples_comp,kw_pair = get_sample_and_kwpair(setting_names=setting_names, smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc)
    fg, ax = plt.subplots(len(kw_pair),len(kw_pair),figsize=(15,15))                
    legend_elements = []
    for i in range(len(samples_comp)):
        param_names = [newProfile(prm_to_cmp)[0] for prm_to_cmp  in kw_pair]
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(setting_name[i],filter=simple_legend)))
        if stnd_lnsprm:
            corner_data = np.array(samples_comp[i]).T 
        elif not ign_lnsprm_qphi:
            corner_dataT = np.array(samples_comp[i])
            ind_e1,ind_e2 = list(kw_pair).index("e1_lens0"),list(kw_pair).index("e2_lens0")
            e1,e2 = corner_dataT[ind_e1],corner_dataT[ind_e2]
            q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
            corner_dataT[ind_e1],corner_dataT[ind_e2] = q,phi
            param_names[ind_e1],param_names[ind_e2] = "$q_1$","$\phi_1$"
            corner_data = corner_dataT.T
        if i==0:
            corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
        else:
            corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    qphi_app = ""
    if not ign_lnsprm_qphi:
        qphi_app = "_qphi"
    if savefig_dir is not None:
        fg.savefig(str(savefig_dir)+"/"+name_pdf+qphi_app+".pdf")
        print("Produced "+name_pdf+qphi_app+".pdf in "+str(savefig_dir))
    else:
        return fg,ax

def get_sample_and_kwpair(setting_names,smpl_mcmc=None,cut_mcmc=None,prm_mcmc=None,backup_path="/backup_results/"):
    if smpl_mcmc is None:
        smpl_mcmc = [get_mcmc_smpl(st,backup_path=backup_path) for st in setting_names]
    if prm_mcmc is None:
        prm_mcmc  = [get_prm_list(st,backup_path=backup_path) for st in setting_names]
    
    all_CP    = all([check_if_CP(st) for st in setting_names])
    all_WWS   = all([not check_if_WS(st) for st in setting_names]) # with source
    param_to_compare= ['theta_E_lens0', 'e1_lens0', 'e2_lens0','theta_E_lens1', 'gamma_ext_lens2','psi_ext_lens2']
    if all_CP:
        param_to_compare = [param_to_compare[0],"gamma_lens0",*param_to_compare[1:] ]

    if all_WWS and not ign_source:
        param_to_compare = [*param_to_compare,'R_sersic_source_light0', 'n_sersic_source_light0']
        
    kw_pair = {}
    for par in param_to_compare:
        pair_i=[]
        for i in range(len(prm_mcmc)):
            if par in prm_mcmc[i]:
                pair_i.append(prm_mcmc[i].index(par))
            else:
                warnings.warn("\nThis parameter is not present in all the settings (maybe fixed?): "+par+"\nIgnored for all settings\n")
                break
        else:
            kw_pair[par] = pair_i
            
    samples_comp = [[] for _ in setting_names]
    for i in kw_pair:
        for j in range(len(smpl_mcmc)):
            smpl_j = np.transpose(smpl_mcmc[j])[kw_pair[i][j]]  # shape: steps
            cut_mcmc_scaled = 0
            if cut_mcmc:
                cut_mcmc_scaled = int(len(smpl_j)*cut_mcmc/1000)
            smpl_j = smpl_j[cut_mcmc_scaled:]
            if "psi_ext" in i:
                samples_comp[j].append(smpl_j*180/np.pi)
            else:
                samples_comp[j].append(smpl_j)               # shape: setting, prm, steps
    prm_comp = list(kw_pair)
    return samples_comp,prm_comp

def plt_SC_LandFP(setting_list,smpl_mcmc=None,prm_mcmc=None,cut_mcmc=None,fermat_mcmc=None,param_fermat=param_names,already_BC=False,stnd_lnsprm=False,
                    col=base_colors,savefig_dir=None,simple_legend=False,name_pdf="lnDf_sup",backup_path="backup_results/"):
    print("Producing superposed corner plot for lens params and Df ")
    settings_names = get_setting_name(setting_list)
    if fermat_mcmc is None:
        fermat_mcmc = [np.array(get_mcmc_fermat(st,backup_path)).tolist() for st in settings_names]
        already_BC  = False
    if not already_BC:
        mcmc_Df = []
        for i in range(len(setting_names)):
            mcmc_DfTi = np.transpose(fermat_mcmc[i])[1:]-np.transpose(fermat_mcmc[i])[0]
            mcmc_Df.append(np.transpose(mcmc_DfTi).tolist())
    else:
        mcmc_Df = fermat_mcmc
    warnings.warn("Given the low S/N of image D, I will discard here and instead consider Delta BC")
    param_fermat[-1] = param_fermat[0].replace("AB","BC")
    param_fermat     = [ p.replace(" [\"^2]","")+" [\"^2]" for p in param_names]
    
    samples_comp,prm_comp  = get_sample_and_kwpair(setting_names=setting_names, smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc)
    _prm_comp = [newProfile(prcmp)[0] for prcmp  in prm_comp] 
    prm_comb2 = [*_prm_comp,*param_fermat]
    
    fg, ax = plt.subplots(len(prm_comb2),len(prm_comb2),figsize=(10,10))
    legend_elements  = []
    for i,sett_nm in enumerate(settings_names):
        col_i = col[i%len(col)]
        if stnd_lnsprm or ign_lnsprm_qphi:
            corner_data = samples_comp[i]
        else:
            corner_data = np.array(samples_comp[i])
            ind_e1,ind_e2 = list(prm_comp).index("e1_lens0"),list(prm_comp).index("e2_lens0")
            e1,e2 = corner_data[ind_e1],corner_data[ind_e2]
            q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
            corner_data[ind_e1],corner_data[ind_e2] = q,phi
            prm_comb2[ind_e1],prm_comb2[ind_e2] = "$q_1$","$\phi_1$"
            corner_data = corner_data.tolist()
        samples_comb2 = np.transpose([*corner_data,*mcmc_Df[i]]).tolist()
        legend_elements.append(Patch(facecolor=col_i,label=strip_setting_name(sett_nm,filter=simple_legend)))
        if i==0:    
            corner(samples_comb2,labels =prm_comb2,show_titles=True,color=col_i,plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
        else:
            corner(samples_comb2,color=col_i,fig=fg,plot_datapoints=False,hist_kwargs= {"density":True})
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
 
    if savefig_dir is not None:
        fg.savefig(str(savefig_dir)+"/"+name_pdf+".pdf")
        print("Produced "+name_pdf+".pdf in "+str(savefig_dir))
    else:
        return fg,ax

        



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot the superposed corner plot of the posterior distribution of the Fermat potential difference and the other lens parameters from the given settings ")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="Cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                        help="Directory name where to save the plot")
    parser.add_argument("-if","--ignore_fermat",dest="ignore_fermat", default=False,action="store_true",
                        help="Ignore superposition of corner plot of fermat potential differences")
    parser.add_argument("-iqp","--ignore_lens_param_qphi",dest="ignore_lens_param_qphi", default=False,action="store_true",
                        help="Ignore superposition of corner plot of lens paramaters (where e1 and e2 have been translate to q and phi)")
    parser.add_argument("-stnd_lnsp","--standard_lens_param",dest="standard_lens_param", default=False,action="store_true",
                        help="Plot the standard superposition of corner plot of lens paramaters (with e1 and e2)")
    parser.add_argument("-sl","--simple_legend",dest="simple_legend", default=False,action="store_true",
                        help="Draw a simplified legend with only the name of the filters")
    parser.add_argument("-is","--ignore_source",dest="ignore_source", default=False,action="store_true",
                        help="Ignore source parameters (if all settings have them)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")

    args = parser.parse_args()
    cut_mcmc        = args.cut_mcmc
    dir_name        = args.dir_name
    setting_names   = args.SETTING_FILES
    ign_Df          = args.ignore_fermat
    stnd_lnsprm     = args.standard_lens_param
    ign_lnsprm_qphi = args.ignore_lens_param_qphi
    simple_legend   = args.simple_legend
    ign_source      = args.ignore_source        

    setting_name    = get_setting_name(setting_names)
    
    #############################
    present_program(sys.argv[0])
    #############################


    backup_path="backup_results"

    save_dir = create_dir_name(setting_names,save_dir="PDF_Sup_Corner",dir_name=dir_name,backup_path=backup_path)
    for sett_name in setting_names:
        svfgp=get_savefigpath(sett_name)
        if svfgp[-1]=="/":
            svfgp = svfgp[:-1]
        os.symlink(os.getcwd()+"/"+svfgp,str(save_dir)+"/.")
    os.system("python check_comment.py "+str(setting_names).replace("[","").replace("]","").replace(",","") +" >> "+str(save_dir)+"/sett_comments.txt &")
    
    smpl_mcmc = [get_mcmc_smpl(st,backup_path=backup_path) for st in setting_names]
    prm_mcmc  = [get_prm_list(st,backup_path=backup_path) for st in setting_names]
    fermat_mcmc = [np.array(get_mcmc_fermat(st,backup_path)).tolist() for st in setting_names]
    already_BC  = False
    if not ign_Df:
        plot_sup_corner_Df(setting_name,savefig_dir=save_dir,fermat_mcmc=fermat_mcmc,simple_legend=simple_legend,param_names=param_names,backup_path=backup_path,already_BC=already_BC) 

    if stnd_lnsprm or not ign_lnsprm_qphi:
        plt_sup_corner_lnsprm(setting_name,smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc,stnd_lnsprm=stnd_lnsprm,ign_lnsprm_qphi=ign_lnsprm_qphi,\
                                savefig_dir=save_dir,simple_legend=simple_legend,backup_path=backup_path)

    plt_SC_LandFP(setting_name,smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,fermat_mcmc=fermat_mcmc,param_fermat=param_names,already_BC=already_BC,stnd_lnsprm=stnd_lnsprm,
                    col=base_colors,savefig_dir=save_dir,simple_legend=simple_legend,backup_path=backup_path)
    
    success(sys.argv[0])


