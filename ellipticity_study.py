#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Test: i want to check the posterior of q,phi of the ellipticity compared to the one obtained from the luminosity of the lens
# instead of e1 e2


# In[3]:


import argparse
import numpy as np
from corner import corner
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from tools import *
from get_res import *
from conversion import qphi_from_e1e2


# In[ ]:





# In[4]:


def get_qphi_post(smpl,prm):
    e1 = smpl.T[prm.index("e1_lens0")]
    e2 = smpl.T[prm.index("e2_lens0")]
    q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
    return q,phi

def get_qphi_post_ll(smpl,prm):
    smpl = get_mcmc_smpl(sets)
    prm = get_mcmc_prm(sets)
    e1 = smpl.T[prm.index("e1_lens_light0")]
    e2 = smpl.T[prm.index("e2_lens_light0")]
    q,phi = qphi_from_e1e2(e1,e2,ret_deg=True)
    return q,phi

def get_ell_prior(sets,npoints=10000):
    sets     = get_setting_module(sets).setting()
    lns_prms_mx = sets.lens_params[-1][0]
    lns_prms_mn = sets.lens_params[-2][0]
    mx_e1    = lns_prms_mx["e1"]
    mx_e2    = lns_prms_mx["e2"]
    mn_e1    = lns_prms_mn["e1"]
    mn_e2    = lns_prms_mn["e2"]    
    e1_prior = np.random.uniform(mn_e1,mx_e1,npoints)
    e2_prior = np.random.uniform(mn_e2,mx_e2,npoints)
    
    return e1_prior,e2_prior


# In[ ]:


if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = argparse.ArgumentParser(description="Produce posterior of q,phi instead of e1,e2 for the main lens")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="setting file(s) to consider")
    args = parser.parse_args()
    
    settings = args.SETTING_FILES
    for sets in settings:
        print("Setting ",get_setting_name(sets))
        svfg_i = get_savefigpath(sets)+"/ellipticity/"
        mkdir(svfg_i)
        smpl = get_mcmc_smpl(sets)
        prm  = get_mcmc_prm(sets)
        q,phi  = get_qphi_post(smpl,prm)
        fg = corner(np.transpose([q,phi]),labels=["q",r"$\phi$"],show_titles=True)
        fg.suptitle("Post. of q,$\phi$")
        plt.tight_layout()
        fg.savefig(svfg_i+"/qphi_post.png")
        print("Created "+svfg_i+"/qphi_post.png")
        if not check_if_SUB(sets):
            q_ll,phi_ll  = get_qphi_post_ll(smpl,prm)
            fg = corner(np.transpose([q_ll,phi_ll]),labels=["q",r"$\phi$"],show_titles=True)
            fg.suptitle("Lens Light")
            plt.tight_layout()
            fg.savefig(svfg_i+"/qphi_post_ll.png")
            print("Created "+svfg_i+"/qphi_post_ll.png")
            plt.close()
            diff_phi = phi - phi_ll
            fg,ax  = plt.subplots()
            ax.hist(diff_phi,histtype="step")
            ax.set_title(r"$\phi_{lns}-\phi_{ll}$")
            plt.tight_layout()
            plt.savefig(svfg_i+"/dphi.png")
            plt.close()
            print("Created "+svfg_i+"/dphi.png")
        plt.close()
        e1_prior,e2_prior = get_ell_prior(sets)
        fg    = corner(np.transpose([e1_prior,e2_prior]),labels=["$e_1^{prior}$",r"$e_2^{prior}$"],show_titles=True)
        fg.savefig(svfg_i+"/e12_prior.png")
        print("Created "+svfg_i+"/e12_prior.png")
        q_prior,phi_prior =  qphi_from_e1e2(e1_prior,e2_prior,ret_deg=True)
        fg = corner(np.transpose([q_prior,phi_prior]),labels=["$q^{prior}$",r"$\phi^{prior}$"],show_titles=True)
        fg.savefig(svfg_i+"/qphi_prior.png")
        plt.close()
        
        
        # compare prior and post for q,phi
        legend_elements  = []
        fg,ax = plt.subplots(2,2,figsize=(6,6))
        truths_qphi = [None,None]
        if check_if_SUB(sets):
            setm = get_setting_module(sets,1)
            truths_qphi = [setm.q_ll,setm.phi_ll]
        corner(np.transpose([q_prior,phi_prior]),truths=truths_qphi,labels=["$q^{prior}$",r"$\phi^{prior}$"],show_titles=True,color="b",hist_kwargs= {"density":True},fig=fg)
        legend_elements.append(Patch(facecolor="b",label="Prior"))
        corner(np.transpose([q,phi]), fig=fg,color="g",hist_kwargs= {"density":True})
        legend_elements.append(Patch(facecolor="g",label="Post. from Lens model"))
        if not check_if_SUB(sets):
            corner(np.transpose([q_ll,phi_ll]),fig=fg,color="r",hist_kwargs= {"density":True})
            legend_elements.append(Patch(facecolor="r",label="Post. from Lens Light model"))
        fg.suptitle(r"Prior of q,$\phi$ vs Post.")
        ax_i = ax[0][1]
        ax_i.legend(handles=legend_elements)
        ax_i.axis("off")
        plt.tight_layout()
        fg.savefig(svfg_i+"/qphi_prior_vs_post.png")
        print("Created "+svfg_i+"/qphi_prior_vs_post.png")
        plt.close()
        
        #TEST -> result from PSO vs mcmc posterior of q,phi
        kw_lens         = get_kwres(sets)["kwargs_results"]["kwargs_lens"]
        e1,e2           = kw_lens[0]["e1"],kw_lens[0]["e2"]
        q_res,phi_res   = qphi_from_e1e2(e1,e2,ret_deg=True)
        print("max(phi),min(phi),phi_res")
        print(max(phi),min(phi),phi_res)
        kw_lens_mcmc    = get_kwres(sets,updated=True)["kwargs_results"]
        e1_mcmc,e2_mcmc = kw_lens_mcmc["e1_lens0"],kw_lens_mcmc["e2_lens0"]
        q_mcmc,phi_mcmc = qphi_from_e1e2(e1_mcmc,e2_mcmc,ret_deg=True)
        print("max(phi),min(phi),phi_mcmc")
        print(max(phi),min(phi),phi_mcmc)
        fg = corner(np.transpose([q,phi]),truths=[q_res,phi_res],labels=["$q^{post}$",r"$\phi^{post}$"],color="b",truth_color="r",show_titles=True,hist_kwargs= {"density":True})
        corner(np.transpose([q,phi]),truths=[q_mcmc,phi_mcmc], color="b",truth_color="g", hist_kwargs= {"density":True},fig=fg)
        fg.suptitle(r"PSO Res.for q,$\phi$ vs Post.")
        plt.tight_layout()
        fg.savefig(svfg_i+"/qphi_post_vs_res.png")
        print("Created "+svfg_i+"/qphi_post_vs_res.png")
                
        #TEST for e1 e2
        fg = corner(np.transpose([e1_prior,e2_prior]),truths=[e1,e2],labels=["$e_1^{prior}$",r"$e_2^{prior}$"],color="b",show_titles=True,hist_kwargs= {"density":True})
        print("e1",max(e1_prior),min(e1_prior),e1)
        print("e2",max(e2_prior),min(e2_prior),e2)
        fg.suptitle(r"PSO Res.for $e_1$,$e_2$ vs Prior")
        plt.tight_layout()
        fg.savefig(svfg_i+"/e12_prior_vs_res.png")
        print("Created "+svfg_i+"/e12_prior_vs_res.png")

        e1_post,e2_post = get_mcmc_smpl_for_prm(sets,"e1_lens0"),get_mcmc_smpl_for_prm(sets,"e2_lens0")

        fg = corner(np.transpose([e1_post,e2_post]),truths=[e1,e2],labels=["$e_1^{post}$",r"$e_2^{post}$"],color="b",truth_color="r",show_titles=True,hist_kwargs= {"density":True})
        corner(np.transpose([e1_post,e2_post]),truths=[e1_mcmc,e2_mcmc], color="b",truth_color="g", hist_kwargs= {"density":True},fig=fg)
        fg.suptitle(r"PSO Res.for $e_1$,$e_2$ vs Post")
        plt.tight_layout()
        fg.savefig(svfg_i+"/e12_post_vs_res.png")
        print("Created "+svfg_i+"/e12_post_vs_res.png")
    success(sys.argv[0])

