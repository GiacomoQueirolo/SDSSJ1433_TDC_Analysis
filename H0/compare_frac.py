# following Schneider and Sluse 2013:
# Source-position transformation: an approximate invariance  in strong gravitational lensing
# section 4.3
"""
Finally, we also display in Fig. 4 the change of the time de-
lay τ/τ̂ normalised by |β|/|β̂|. While the delays should scale like
the source positions in the case of a MST, they do not do so
for the SPT. The ratio (τ/τ̂)/(|β|/|β̂|) is almost constant but dif-
fers significantly from 1. In addition, we see that for sources
interior to the inner (astroid) caustic, the time delay ratios be-
tween lensed images are not conserved. The error bars in Fig. 4
show, for each source, the spread of τ/τ̂ among the lensed im-
ages. This spread may reach 10% for some peculiar source po-
sitions, but generally arises from the deviation of only one of
the lensed images. Owing to the accuracies on the time delays
of a few percents currently achieved for quads (e.g. Fassnacht
et al. 2002; Courbin et al. 2011; Tewes et al. 2013), time delay
measurements in quadruply lensed systems could in the most
favourable cases break the SPT independently of H0 .
"""
# I want to compare the fraction of dtai/dtaj with dphi_ai/dphi_aj with the rel.
# error, and see if they are compatible

import os,sys
import numpy as np
import argparse as ap

from Utils.tools import *
from Utils.statistical_tools import marginalise_prob

from Plots.plotting_tools import get_median_with_error
import matplotlib.pyplot as plt
from H0.H0_Combined_reworked import kwpycs,kwlnst
from H0.H0_Data import get_Dt_post,get_Df_post,get_kwdt,verify_BC,create_PH0resdir

def _get_ratio(dx,sig_dx):
    rt_dx   = []
    sig_rdx = []
    for i in range(len(dx)): #ab,ac,ad
        rt_dx_i=[]
        sg_dx_i=[]
        for j in range(len(dx)):
            if i!=j:
                dxij = dx[i]/dx[j]
                rt_dx_i.append(dxij)
                sg_dx_i.append(dxij*np.sqrt(np.sum([(sig_dx[k]/dx[k])**2 for k in [i,j]])))
        rt_dx.append(rt_dx_i)
        sig_rdx.append(sg_dx_i)

    return rt_dx,sig_rdx

def _get_ratio_labels(labels):
    ratio_labels = []
    for i in range(len(labels)): #ab,ac,ad
        ratio_labels_i= []
        for j in range(len(labels)):
            if i!=j:
                ratio_labels_i.append(str(labels[i])+"/"+str(labels[j]))
        ratio_labels.append(ratio_labels_i)
    print(ratio_labels)
    return ratio_labels

def get_ratios(df,dt,sig_df,sig_dt,labels):
    assert len(dt)==len(df) # must be same length
    rt_df,sig_rdf = _get_ratio(df,sig_df)
    rt_dt,sig_rdt = _get_ratio(dt,sig_dt)
    ratio_labels  = _get_ratio_labels(labels=labels)

    return rt_df,sig_rdf,rt_dt,sig_rdt,ratio_labels

def plot_ratios(rdt,sig_rdt,rdf,sig_rdf,ratio_labels,path="."):
    assert len(rdt)==len(rdf) # must be same length
    
    fig,ax = plt.subplots(len(rdt),len(rdt[1]))
    ind = 1
    for i in range(len(rdt)):
        for j in range(len(rdt[i])):
            ax[i][j].scatter(rdf[i][j],ind,label="Df_"+ratio_labels[i][j],color="r")
            ax[i][j].errorbar(rdf[i][j],ind,xerr=sig_rdf[i][j],fmt=".r")
            ax[i][j].scatter(rdt[i][j],ind+0.3,label="Dt_"+ratio_labels[i][j],color="b")
            ax[i][j].errorbar(rdt[i][j],ind+0.3,xerr=sig_rdt[i][j],fmt=".b")
            ax[i][j].set_title(ratio_labels[i][j])
            ax[i][j].legend()

    figname = path+"/Ratio_Compared.pdf"
    print("Saving "+figname)
    fig.tight_layout()
    fig.savefig(figname)
    return 0

if __name__ == '__main__':
        
    present_program(sys.argv[0])
    
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Combine the posterior for time delay and Fermat potential differences at the images position to constrain the Hubble parameter H0",
                               formatter_class=ap.RawTextHelpFormatter)
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"
    #help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of savedir for multiplied posteriors)"
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of combined setting file)"
    help_plotdir   = "Name of the directory which will be containing the plot(should be same as post_dir)"
    help_overwrite = "Overwrite previous result. If False and I find the same result directory, I will create a new one with today's date."
    help_margOm    = "Consider marginalised Posterior with respect to Omega_m"
    parser.add_argument("-dplot","--dir_plot",dest='dir_plot', type=str,default="PH0",
                        metavar='dir_plot', action='store',
                        help=help_plotdir)
    parser.add_argument("-dtn","--Dt_name",dest='dt_name', type=str,default=".",
                        metavar='dt_name', action='store',
                        help=help_timedelay)
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default=".",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
    args      = parser.parse_args()
    dir_plot  = args.dir_plot

    dt_name   = args.dt_name
    df_name   = args.df_name
    # loading data
    conf_dt,Dt_res = get_Dt_post(dt_name=dt_name,link=False,**kwpycs)
    combined_setting,PDF_Df,PDF_Df_bins = get_Df_post(df_name=df_name,link=False,**kwlnst)
    verify_BC(conf_dt=conf_dt,comb_sett=combined_setting)
    kwargs_dt = get_kwdt(Dt_res)    




    bin_densities = marginalise_prob(PDF_Df,PDF_Df_bins) 
    res_df,sig_df = [],[]

    for i in range(3):
        for j in range(3):
            if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                if i==j:
                    results_Dfi = get_median_with_error(bin_densities[i],PDF_Df_bins[i],ret_str=False)
                    res_df.append(results_Dfi[0])
                    sig_df.append(np.mean(results_Dfi[1]))
    sig_dt = np.diag(kwargs_dt["cov"])**2
    lbls = [ lbl.replace(r"$\Delta\phi_{","").replace(r"}$","") for lbl in combined_setting.labels_df]
    rt_df,sig_rdf,rt_dt,sig_rdt,ratio_labels = get_ratios(df=res_df,dt=kwargs_dt["mean"],sig_df=sig_df,sig_dt=sig_dt,labels=lbls)

    res_dir   = "./results/"
    PH0_resdir = create_PH0resdir(res_dir=res_dir,dir_ph0=dir_plot,verbose=True,overwrite=False)

    plot_ratios(rt_df,sig_rdf,rt_dt,sig_rdt,ratio_labels,path=PH0_resdir)
    success(sys.argv[0])
