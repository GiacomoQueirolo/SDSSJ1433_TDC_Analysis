import os,sys
import corner
import argparse
import warnings
import numpy as np
from corner import corner
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from tools import *
from get_res import *
from Prior import get_prior
from conversion import qphi_from_e1e2
from plotting_tools import base_colors
from get_param_title import newProfile
from source_pos import get_source_gen,get_source_pos_PSO

def get_prior_src_pos(sett,npoints=1000,saveprior=False):
    # this function is not used here (yet) but might be useful
    # from priors get mcmc prior of source pos
    sett = get_setting_module(sett).setting()
    mcmc_prior,param_mcmc = get_prior(sett,npoints)
    MCP_src_ra,MCP_src_dec = [],[] #MCMC prior source ra and dec
    for i in range(len(mcmc_prior)):
        kwres_i = conv_mcmc_i_to_kwargs(sett,mcmc_prior[i])
        kw_lens = kwres_i["kwargs_lens"]
        ra_im   = kwres_i["kwargs_ps"][0]["ra_image"]
        dec_im  = kwres_i["kwargs_ps"][0]["dec_image"]
        ra_src,dec_src = get_source_gen(ra_im,dec_im,kw_lens,sett)
        MCP_src_ra.append(*ra_src)
        MCP_src_dec.append(*dec_src)
    if saveprior:
        savepath = get_savemcmcpath(sett)
        save_json([MCP_src_ra,MCP_src_dec],savepath)
    return MCP_src_ra,MCP_src_dec
    
    

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Simply plot the posterior of the source wrt the other params")
    parser.add_argument("-PP", "--plot_prior", action="store_true", dest="plot_prior", default=False,
                        help="Plot the prior distribution for the source parameters")
    parser.add_argument("-np", "--npoints", type=int,dest="npoints", default=1000,
                        help="Number of points for the MCMC prior")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=".",
                    help="Directory name where to save the plot")
    parser.add_argument("-sl","--simple_legend",dest="simple_legend", default=False,action="store_true",
                    help="Draw a simplified legend with only the name of the filters")
    parser.add_argument("-ilp","--ignore_lens_params",dest="ignore_lens_params", default=False,action="store_true",
                    help="Ignore superposition of corner plot of lens paramaters (only source params)")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args          = parser.parse_args()
    npoints       = args.npoints
    plot_prior    = args.plot_prior
    settings      = args.SETTING_FILES
    dir_name      = args.dir_name
    simple_legend = args.simple_legend
    ign_lsn_prms  = args.ignore_lens_params

    backup_path      = "backup_results"
    prior_name       = "mcmc_prior.png"
    smpl_mcmc_prior  = []
    prm_mcmc_prior   = []
    for sets in settings:
        savefig_path = get_savefigpath(sets)
        print_setting(sets)
        mcmc_prior,param_mcmc = get_prior(sets,npoints)
        if plot_prior:
            savefig_path = get_savefigpath(sets,backup_path)
            lbls = [newProfile(prm)[0] for prm in param_mcmc]
            corner(mcmc_prior,labels=lbls,show_titles=True)
            plt.savefig(savefig_path+"/"+prior_name)
            plt.close()
            src_pr_x,src_pr_y = get_prior_src_pos(sett,npoints=1000,saveprior=False)
            kw_source_pso = get_source_pos_PSO(sets)
            src_pso_x,src_pso_y = kw_source_pso["source_ra"],kw_source_pso["source_dec"]
            corner(np.array([src_pr_x,src_pr_y],lables=["Ra prior Source","Dec prior Source"],truths=[src_pso_x,src_pso_y]))
            plt.savefig(savefig_path+"/Prior_source_pos.png")
        smpl_mcmc_prior.append(mcmc_prior)
        prm_mcmc_prior.append(param_mcmc)
        
    if len(settings)>1:
        save_dir_sup = create_dir_name(settings,save_dir="PDF_Sup_Corner/priors/",dir_name=dir_name,backup_path=backup_path,copy_settings=True)
        
        col = base_colors

        all_CP    = all([check_if_CP(st) for st in settings])
        all_WWS   = all([not check_if_WS(st) for st in settings]) # with source
        param_to_compare = []
        if not ign_lsn_prms:
            param_to_compare = ['theta_E_lens0', 'e1_lens0', 'e2_lens0','theta_E_lens1', 'gamma_ext_lens2','psi_ext_lens2']
            if all_CP:
                param_to_compare = [param_to_compare[0],"gamma_lens0",*param_to_compare[1:] ]
        
        if all_WWS:
            param_to_compare = [*param_to_compare,'R_sersic_source_light0', 'n_sersic_source_light0', 'e1_source_light0', 'e2_source_light0']
        
        kw_pair = {}
        for par in param_to_compare:
            pair_i=[]
            for i in range(len(prm_mcmc_prior)):
                if par in prm_mcmc_prior[i]:
                    pair_i.append(prm_mcmc_prior[i].index(par))
                else:
                    warnings.warn("\nThis parameter is not present in all the settings (maybe fixed?): "+par+"\nIgnored for all settings\n")
                    break
            else:
                kw_pair[par] = pair_i
            
        samples_comp = [[] for i in settings]
        for i in kw_pair:
            for j in range(len(smpl_mcmc_prior)):
                smpl_j = np.transpose(smpl_mcmc_prior[j])[kw_pair[i][j]]  # shape: steps
                if "psi_ext" in i:
                    samples_comp[j].append(smpl_j*180/np.pi)
                else:
                    samples_comp[j].append(smpl_j)               # shape: setting, prm, steps
        
        fg, ax = plt.subplots(len(kw_pair),len(kw_pair),figsize=(15,15))                
        legend_elements = []
        for i in range(len(samples_comp)):
            param_names = [newProfile(prm_to_cmp)[0] for prm_to_cmp  in kw_pair]
            legend_elements.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(settings[i],filter=simple_legend)))
            #if stnd_lnsprm:
            #    corner_data = np.array(samples_comp[i]).T 
            #elif not ign_lnsprm_qphi:
            corner_dataT = np.array(samples_comp[i])
            if not ign_lsn_prms:
                ind_e1,ind_e2 = list(kw_pair).index("e1_lens0"),list(kw_pair).index("e2_lens0")
                e1,e2 = corner_dataT[ind_e1],corner_dataT[ind_e2]
                q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
                corner_dataT[ind_e1],corner_dataT[ind_e2] = q,phi
                param_names[ind_e1],param_names[ind_e2] = "$q_1$","$\phi_1$"
            corner_data = corner_dataT.T
            if i==0:
                corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg)
            else:
                corner(corner_data,color=col[i%len(col)],fig=fg,plot_datapoints=False,hist_kwargs= {"density":True})
        axdel=ax[0][-1]
        axdel.legend(handles=legend_elements)
        axdel.axis("off")
        fg.savefig(str(save_dir_sup)+"/sup_prior.pdf")
        print("Produced sup_prior.pdf in "+str(save_dir_sup))


    success(sys.argv[0])
