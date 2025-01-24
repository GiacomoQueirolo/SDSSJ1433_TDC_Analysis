import pickle
import os, sys,gc
import argparse
import warnings
import numpy as np
from corner import corner
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

#my libs
from Utils.tools import *
from Utils.get_res import *
from Utils.get_param_title import newProfile_w_images as newProfile

from Data.Param import get_prm_list
from Data.conversion import qphi_from_e1e2

from Plots.plotting_tools import base_colors
from Posterior_analysis.fermat_pot_analysis import get_mcmc_Df
#from Posterior_analysis.fermat_pot_analysis import labels_Df as  param_names   
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
"""
def plot_sup_corner_Df_old(setting_list,fermat_mcmc=None,col=base_colors,savefig_dir=None,name_pdf="Df_sup",cut_mcmc=None,simple_legend=False,
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
""" 
    def plot_sup_corner_Df_1xtime(setting_list,Prev_Fig=None,param_names=None,BC=True,cut_mcmc=None,
                       col=base_colors,savefig_dir=None,name_pdf="Df_sup",simple_legend=False,
                       backup_path="backup_results/",_produce_plot=True):
    print("Producing Df superposed corner plot, one at the time")
    #print("Note: imput is assumed to be already the difference of fermat Potentials") #doesn't matter as we also get the param_names
    # for fermat potentials
    settings_name = get_setting_name(setting_list)

    if Prev_Fig is None:
        Prev_Fig, ax = plt.subplots(3,3,figsize=(10,10))

        if param_names is None:
            if BC:
                from Posterior_analysis.fermat_pot_analysis import labels_Df_BC as param_names   
            else:
                from Posterior_analysis.fermat_pot_analysis import labels_Df_AD as param_names   
            param_names = [ p.replace(r" [arcsec$^2$]","")+r" [arcsec$^2$]" for p in param_names]
        lgnd_hndls = []
    else:
        lgnd_hndls = Prev_Fig.axes[2].legend_.legendHandles
        
    for i,sett in enumerate(settings_name):
        try:
            Df_mcmc_i = get_mcmc_Df(sett,backup_path,noD=BC)
            mcmc_Dfi  = cut_mcmc(Df_mcmc_i,cut=cut_mcmc,scale=1000)

            lgnd_hndls.append(Patch(facecolor=col[i%len(col)],label=strip_setting_name(settings_name[i],filter=simple_legend)))
            Prev_Fig = corner(mcmc_Dfi,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=Prev_Fig,title_fmt=None)
            Prev_Fig.axes[2].legend(handles=[*lgnd_hndls])
            Prev_Fig.axes[2].axis("off")
            
            name_corner_obj = "_".join([lgd.get_label() for lgd in lgnd_hndls])+".fig.pkl"
            with open(f"{savefig_dir}/{name_corner_obj}","wb") as f:
                pickle.dump(Prev_Fig,f)

        except MemoryError:
            print("Memory error, going iteratively")
            del Df_mcmc_i
            del mcmc_Dfi
            gc.collect()

            Prev_Fig = plot_sup_corner_Df(setting_list=setting_list[i-1:],Prev_Fig=Prev_Fig,param_names=param_names,
                         BC=BC,cut_mcmc=cut_mcmc,col=col,savefig_dir=savefig_dir,name_pdf=name_pdf,
                         simple_legend=simple_legend,backup_path=backup_path,_produce_plot=False)

    if savefig_dir and _produce_plot:
        Prev_Fig.savefig(str(savefig_dir)+"/"+name_pdf+".pdf")
        print("Produced "+name_pdf+".pdf in "+str(savefig_dir))
    else:
        return Prev_Fig


def plot_sup_corner_Df(setting_list,Df_mcmc=None,param_names=None,BC=True,cut_mcmc=None,
                        col=base_colors,savefig_dir=None,name_pdf="Df_sup",simple_legend=False,
                        backup_path="backup_results/",Prev_Fig=None):
    print("Producing Df superposed corner plot")
    #print("Note: imput is assumed to be already the difference of fermat Potentials") #doesn't matter as we also get the param_names
    # for fermat potentials
    settings_name = get_setting_name(setting_list)
    if Df_mcmc is None:
        Df_mcmc = [get_mcmc_Df(st,backup_path,noD=BC) for st in settings_name]
    if param_names is None:
        if BC:
            warnings.warn("Given the low S/N of image D, I will discard here and instead consider Delta BC")
            from Posterior_analysis.fermat_pot_analysis import labels_Df_BC as param_names   
        else:
            from Posterior_analysis.fermat_pot_analysis import labels_Df_AD as param_names   
        param_names = [ p.replace(r" [arcsec$^2$]","")+r" [arcsec$^2$]" for p in param_names]
    if Prev_Fig:
        fg = Prev_Fig
        ax = fg.axes
    else:
        fg, ax = plt.subplots(3,3,figsize=(10,10))
    legend_elements  = []
    for i,mcmc_Dfi in enumerate(Df_mcmc):
        cut_mcmc_scaled = 0
        mcmc_DfiT =  np.transpose(mcmc_Dfi)
        if cut_mcmc:
            cut_mcmc_scaled = int(len(mcmc_DfiT)*cut_mcmc/1000)
        mcmc_Dfi   = np.array(mcmc_Dfi[cut_mcmc_scaled:])
        lbl = str(strip_setting_name(settings_name[i],filter=simple_legend)).upper()
        if simple_legend:
            lbl = r"$"+lbl.replace("_","_{")+r"}$"
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=lbl))
        if i==0 :
            corner(mcmc_Dfi,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg,title_fmt=None)
        else:
            corner(mcmc_Dfi,color=col[i%len(col)],fig=fg,plot_datapoints=False,hist_kwargs= {"density":True},title_fmt=None)
    axdel=ax[0][2]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
    #plt.tight_layout()
    del mcmc_Dfi
    if savefig_dir:
        fg.savefig(str(savefig_dir)+"/"+name_pdf+".pdf")
        print("Produced "+name_pdf+".pdf in "+str(savefig_dir))
    else:
        return fg,ax

def plt_sup_corner_lnsprm(setting_list,smpl_mcmc=None,prm_mcmc=None,cut_mcmc=None,stnd_lnsprm=False,ign_lnsprm_qphi=False,
                    col=base_colors,savefig_dir=None,simple_legend=False,name_pdf="lnsprm_sup",backup_path="backup_results/",radec=False):
    print("Producing lens params superposed corner plot")
    setting_name = get_setting_name(setting_list) 
    samples_comp,kw_pair = get_sample_and_kwpair(setting_names=setting_name, smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc,radec=radec)
    fg, ax = plt.subplots(len(kw_pair),len(kw_pair),figsize=(2*len(kw_pair),2*len(kw_pair)))                
    legend_elements = []
    for i in range(len(samples_comp)):
        param_names = [newProfile(prm_to_cmp)[0] for prm_to_cmp  in kw_pair]
        lbl = strip_setting_name(setting_name[i],filter=simple_legend)
        if simple_legend:
            lbl = r"$"+strip_setting_name(setting_name[i],filter=simple_legend).replace("_","_{").upper()+r"}$"
        legend_elements.append(Patch(facecolor=col[i%len(col)],label=lbl))
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
            """if test:
                ind_ra1,ind_dec1 = list(kw_pair).index("ra"),list(kw_pair).index("dec")
                ra,dec = corner_dataT[ind_ra1:ind_ra1+4],corner_dataT[ind_dec1:ind_dec1+4]
                ra_m,dec_m  = np.mean(ra,axis=0),np.mean(dec,axis=0)
                img_order = image_order(ra_m,dec_m,ret_order=True)
                nm_img = "A","B","C","D"
                ord_nm_img = nm_img[img_order]
                for i in range(ind_ra1,ind_ra1+4):
                    param_names[i] = r"$ra_"+ord_nm_img[i]+"$"
                for i in range(ind_dec1,ind_dec1+4):
                    param_names[i] = r"$dec_"+ord_nm_img[i]+"$"
            """
        if i==0:
            corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg,title_fmt=None)
        else:
            corner(corner_data,labels=param_names,show_titles=True,color=col[i%len(col)],plot_datapoints=False,hist_kwargs= {"density":True},fig=fg,title_fmt=None)
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

def get_sample_and_kwpair(setting_names,smpl_mcmc=None,cut_mcmc=None,prm_mcmc=None,backup_path="/backup_results/",radec=False):
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
        
    if radec:
        param_to_compare = [*param_to_compare,"center_x_lens0","center_y_lens0","cnt_pert_x","cnt_pert_y"]
        param_to_compare = [*param_to_compare,*["ra_image_"+str(i) for i in range(4)],*["dec_image_"+str(i) for i in range(4)]]

    kw_pair = {}
 
    ra_i,dec_i = 0,0
    for par in param_to_compare:
        pair_i=[]
        for i in range(len(prm_mcmc)):
            if par in prm_mcmc[i]:
                pair_i.append(prm_mcmc[i].index(par))
            elif radec and "image" in par:
                if "ra" in par and "ra_image" in prm_mcmc[i]:
                    first_im = prm_mcmc[i].index("ra_image")
                    im_indx = first_im+int(par[-1])#ra_i
                    ra_i+=1
                elif "dec" in par and "dec_image" in prm_mcmc[i]:
                    first_im = prm_mcmc[i].index("dec_image")
                    im_indx = first_im+int(par[-1])
                    dec_i+=1
                pair_i.append(im_indx)
            elif radec and "cnt_pert" in par:
                if "center_x_lens_light1" in prm_mcmc[i]:
                    if "x" in par:
                        par = "center_x_lens1"
                        pair_i.append(prm_mcmc[i].index("center_x_lens_light1"))
                    elif "y" in par:
                        par = "center_y_lens1"
                        pair_i.append(prm_mcmc[i].index("center_y_lens_light1"))
                    else:
                        raise RuntimeError("something wrong A")
                elif "center_x_lens_light0" in prm_mcmc[i] and not "center_x_lens_light1" in prm_mcmc[i]:
                    if "x" in par:
                        par = "center_x_lens1"
                        pair_i.append(prm_mcmc[i].index("center_x_lens_light0"))
                    elif "y" in par:
                        par = "center_y_lens1"
                        pair_i.append(prm_mcmc[i].index("center_y_lens_light0"))
                    else:
                        raise RuntimeError("something wrong B")
                else:
                    raise RuntimeError("something wrong C")
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

def plt_SC_LandFP(setting_list,smpl_mcmc=None,prm_mcmc=None,cut_mcmc=None,Df_mcmc=None,param_fermat=None,BC=True,stnd_lnsprm=False,
                    col=base_colors,savefig_dir=None,simple_legend=False,name_pdf="lnDf_sup",backup_path="backup_results/",radec=True):
    print("Producing superposed corner plot for lens params and Df ")
    settings_names = get_setting_name(setting_list)
    if Df_mcmc is None:
        Df_mcmc = [np.transpose(get_mcmc_Df(st,backup_path,noD=BC)) for st in setting_list]
    if param_fermat is None:
        if BC:
            from Posterior_analysis.fermat_pot_analysis import labels_Df_BC as param_fermat   
        else:
            from Posterior_analysis.fermat_pot_analysis import labels_Df_AD as param_fermat   
        param_fermat = [ p.replace(r" [arcsec$^2$]","")+r" [arcsec$^2$]" for p in param_fermat]
            
    samples_comp,prm_comp  = get_sample_and_kwpair(setting_names=settings_names, smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc,radec=radec)
    _prm_comp = [newProfile(prcmp)[0] for prcmp  in prm_comp] 
    prm_comb2 = [*_prm_comp,*param_fermat]
    
    fg, ax = plt.subplots(len(prm_comb2),len(prm_comb2),figsize=(len(prm_comb2)*3,len(prm_comb2)*3))
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
        samples_comb2 = np.transpose([*corner_data,*Df_mcmc[i]])
        lbl = strip_setting_name(sett_nm,filter=simple_legend)
        if simple_legend:
            lbl = r"$"+strip_setting_name(sett_nm,filter=simple_legend).replace("_","_{").upper()+r"}$"
        legend_elements.append(Patch(facecolor=col_i,label=lbl))
        if i==0:    
            corner(samples_comb2,labels =prm_comb2,show_titles=True,color=col_i,plot_datapoints=False,hist_kwargs= {"density":True},fig=fg,title_fmt=None)
        else:
            corner(samples_comb2,color=col_i,fig=fg,plot_datapoints=False,hist_kwargs= {"density":True},title_fmt=None)
    axdel=ax[0][-1]
    axdel.legend(handles=legend_elements)
    axdel.axis("off")
 
    if savefig_dir is not None:
        fg.savefig(str(savefig_dir)+"/"+name_pdf+".pdf")
        print("Produced "+name_pdf+".pdf in "+str(savefig_dir))
    else:
        return fg,ax

"""
def plt_SC_LandFP_old(setting_list,smpl_mcmc=None,prm_mcmc=None,cut_mcmc=None,fermat_mcmc=None,param_fermat=param_names,already_BC=False,stnd_lnsprm=False,
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

"""



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Plot the superposed corner plot of the posterior distribution of the Fermat potential difference and the other lens parameters from the given settings ")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="Cut the first <c> %% of steps of the mcmc to ignore them")
    parser.add_argument("-n","--name",type=str,dest="dir_name", default=None,
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
    parser.add_argument("-AD", dest="AD", default=False,action="store_true",
                        help="Consider AD couple instead of BC")
    parser.add_argument("-rd","--radec", dest="radec", default=False,action="store_true",
                        help="Also plot ra-dec image coordinates")
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
    BC              = not args.AD
    radec           = args.radec
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
        if not os.path.lexists(str(save_dir)+"/."):
            print("linking: ln -s"+os.getcwd()+"/"+svfgp, str(save_dir)+"/.")
            os.symlink(os.getcwd()+"/"+svfgp, str(save_dir)+"/.")
    os.system("python check_comment.py "+str(setting_names).replace("[","").replace("]","").replace(",","") +" >> "+str(save_dir)+"/sett_comments.txt &")
    
    smpl_mcmc = None#[get_mcmc_smpl(st,backup_path=backup_path) for st in setting_names]
    prm_mcmc  = None#[get_prm_list(st,backup_path=backup_path) for st in setting_names]
    Df_mcmc   = None#[np.array(get_mcmc_Df(st,backup_path,noD=BC)).tolist() for st in setting_names]
    
    if not ign_Df:
        """try:
            Df_mcmc = [get_mcmc_Df(st,backup_path=backup_path) for st in setting_names]
        except MemoryError:
            nsteps = 100000
            print(f"MemoryError: Taking only the last {nsteps} of the chain:")
            Df_mcmc = []
            for sett in get_setting_module(setting_names):
                mcmc_path = get_savemcmcpath(sett)+"/mcmc_Df_"+get_setting_name(sett).replace("settings_","").replace(".py","")+".json"
            with open(mcmc_path,"r") as f:
                Df_mcmc.append(json.load(f)[-nsteps:])"""
        try:
            plot_sup_corner_Df(setting_name,savefig_dir=save_dir,Df_mcmc=Df_mcmc,BC=BC,simple_legend=simple_legend,backup_path=backup_path)
        except MemoryError:
            print("Memory error, using Df plot_sup_corner_Df_1xtime")
            plot_sup_corner_Df_1xtime(setting_name,savefig_dir=save_dir,BC=BC,simple_legend=simple_legend,backup_path=backup_path)
    gc.collect()
    try:
        smpl_mcmc = [get_mcmc_smpl(st,backup_path=backup_path) for st in setting_names]
    except MemoryError:
        nsteps = 100000
        print(f"MemoryError: Taking only the last {nsteps} of the chain:")
        smpl_mcmc = []
        for sett in get_setting_module(setting_names):
            mcmc_path = get_savemcmcpath(sett)+"/mcmc_smpl_"+get_setting_name(sett).replace("settings_","").replace(".py","")+".json"
            with open(mcmc_path,"r") as f:
                smpl_mcmc.append(json.load(f)[-nsteps:])
    if stnd_lnsprm or not ign_lnsprm_qphi:
        plt_sup_corner_lnsprm(setting_name,smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,cut_mcmc=cut_mcmc,stnd_lnsprm=stnd_lnsprm,ign_lnsprm_qphi=ign_lnsprm_qphi,\
                                savefig_dir=save_dir,simple_legend=simple_legend,backup_path=backup_path,radec=radec)

    if not ign_Df:
        plt_SC_LandFP(setting_name,smpl_mcmc=smpl_mcmc,prm_mcmc=prm_mcmc,Df_mcmc=Df_mcmc,BC=BC,stnd_lnsprm=stnd_lnsprm,
                    col=base_colors,savefig_dir=save_dir,simple_legend=simple_legend,backup_path=backup_path,radec=radec)
    
    success(sys.argv[0])



