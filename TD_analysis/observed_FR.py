#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

# get the observed Flux Ratio
import sys,os
import pickle
import importlib
import numpy as np
import argparse as ap
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from Utils.tools import *

"""
# soon the function FR_kw should be ignored
# for a more precise analysis
"""
def FR_kw(config,analysis_path,savefig_path):
    all_FR = []
    naming_FR = []
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = analysis_path/str("analysis_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            if mltype_i=="splml":
                continue
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                mlfp_str="_mlfp_"+str(mlfp)
                naming = str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)
                saveml_path = saveknt_path/naming
                with open(str(saveml_path)+"/tot_mag_shift.data","rb") as f:
                    mag_shift= pickle.load(f)
                mag_shifts = np.array(mag_shift[1:])- mag_shift[0] #rel to A 
                # mag_shifts are the shift TO BE APPLIED to the lcs in order to superpose them
                # as in lcs_1 + mag_shifts[1] = lcs_2 + mag_shifts[2]
                # THEREFORE the magnitude shift due to the lens is the opposite of that
                # that is mag_lens = - mag_shifts
                # and since F/F' = 10^(-(m-m')/2.5) -> FR = F/F_A(due to lensing) = 10^(-mag_lens/2.5) =\
                # FR = 10^(mag_shifts/2.5) 
                FR = 10**(mag_shifts/2.5) #F_i/F_A
                all_FR.append(FR)
                naming_FR.append(naming)
    all_FR = np.transpose(all_FR).tolist()
    FR_med = np.median(all_FR,axis=1)

    kwargs_FR = {"FR":all_FR,"FR_names":naming_FR,"FR_median":FR_med}
    with open(savefig_path/str("FR_res.pkl"),"wb") as f:
        pickle.dump(kwargs_FR,f)
    return kwargs_FR


# In[ ]:


def plot_FR_kw(config,savefig_path,kwargs_FR):
    FR = kwargs_FR["FR"]
    FR_med = kwargs_FR["FR_median"]
    #FR_names=kwargs_FR["FR_names"]
    fg,axes= plt.subplots(1,3,figsize=(15,15)) 
    fntsz = 20

    colors= config.colors

    def fmt_new(y,pos):
        return ""

    FR_names = ["AB","AC"]

    for j,ax in enumerate(axes[:2]):
        ax.yaxis.set_major_formatter(FuncFormatter(fmt_new))
        if j==1:
            ax.set_title("Flux Ratio obtained from lcs analysis",fontsize=fntsz)

        y_rnd = [0.5 + j for j in range(len(all_FR[j])+1,0,-1)]

        for i,fri in enumerate(all_FR[j]):
            colj=colors[int(i%len(colors))]
            ax.scatter(fri ,y_rnd[i], c=colj,label=naming_FR[j])
            #ax.text(fri-.05*len(str(np.round(fri,3))),y_rnd[i]+.2,s=str(np.round(fri,2)),c=colj,fontsize=fntsz)
            ax.text(fri,y_rnd[i],s=str(np.round(fri,2)),c=colj,fontsize=fntsz)

        ax.set_xlabel(r"$FR_{"+FR_names[j]+"}$ []",fontsize=fntsz)
        ax.set_ylim(min(y_rnd)-1,max(y_rnd)+1)

        ax.scatter(FR_med[j] ,y_rnd[-1], c="k",label="Median")
        #ax.text(FR_med[j]-.05*len(str(np.round(FR_med[j],3))),y_rnd[-1]+.2,s=str(np.round(FR_med[j],2)),c="k",fontsize=fntsz)
        ax.text(FR_med[j],y_rnd[-1],s=str(np.round(FR_med[j],2)),c="k",fontsize=fntsz)

    ax_del = axes[2]
    ax_del.axis("off")
    ln_lgnd = []
    for j in range(len(all_FR[-1])):
        if j >= len(colors):
            colj=colors[int(j%len(colors))]
        else:
            colj=colors[j]
        y_lbl = naming_FR[j] 
        ln, = ax_del.plot(1,1,c=colj,label=y_lbl)
        ln_lgnd.append(ln)


    ln, = ax_del.plot(1,1,c="k",label="Median")
    ln_lgnd.append(ln)

    #if len(ln_lgnd)<20:
    ax_del.legend(fontsize=17)   
    #else:
    #    set_lgnd = ax_del.legend(handles=ln_lgnd[:int(len(ln_lgnd)/2.)],loc="upper left")
    #    ax_del.add_artist(set_lgnd)
    #    ax_del.legend(handles=ln_lgnd[int(len(ln_lgnd)/2.):],loc="upper center")



    plt.tight_layout()
    #plt.legend(fontsize=fntsz)
    plt.savefig(savefig_path/str("FR_obs.png"))


# In[ ]:


################WIP#############################


# In[ ]:


# 20th Sep 22
# consider the same error propagation of dt but for dmag


# load the run result
# diff of dmag from truemagshift and magshist -> problem
from TD_analysis import stnd_plot
from TD_analysis.stnd_handling_data import *

def get_res_Group_mag_i(config,knt,mltype,ml_config):
    """
    Get the result Group
    """
    
    savefig_path = pth.Path(config.analysis_directory)
    saveknt_path = savefig_path/str("analysis_kn"+str(knt))
            
    mllist_name,mlfp = ml_config 
    mlfp_str  = "_mlfp_"+str(mlfp)
    data_path = str(saveknt_path/str("ml"+mltype[:-2]+mllist_name+mlfp_str) )
    sim_path  = data_path.replace("Analysis","Simulation").replace("analysis","simulation")
    sim_path  = sim_path+"/sims_" + config.simset_mock+"_opt_"+config.optset[0]

    #ERROR
    error = Error_mag(sim_path)
    error.create_error()
    #RES
    res_Group = getresults_mag(data_path,error=error,labels=config.delay_labels)
    return res_Group

def combine_models_mag(config,sigmathresh=None):
    if sigmathresh is None:
        sigmathresh = config.sigmathresh
    wddir = config.combined_directory 
    marginalisation_plot_dir = wddir+'/figure_MAG/marginalisation_plots/'
    mkdir(marginalisation_plot_dir)        

    indiv_marg_dir = marginalisation_plot_dir + config.name_marg_spline + '/'
    mkdir(indiv_marg_dir)        

    marginalisation_dir = wddir+ config.name_marg_spline + '/'
    mkdir(marginalisation_dir)        

    colors = config.colors
    color_id = 0

    group_list = []

    opt = config.optset[0]
        
    # Redefine the keyword here because you don't necessary want to marginalise over everything
    # mlknotsteps_marg now correspond to ml_config, which therefore contains all combinations of ml_lc and degrees

    masked_A = config.maskA
    
    savefig_path=pth.Path(config.analysis_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            if mltype_i!="splml":
                
                for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                    mllist_name,mlfp = ml_config 
                    group_i = get_res_Group_mag_i(config,knt,mltype_i,ml_config)
                    group_i.color=colors[color_id]
                    name = config.combkw[knt_i][mlt_i][mlc_i].replace("_"," ")
                    group_i.name = name
                    group_list.append(group_i)
                    color_id += 1
                    if color_id >= len(colors):
                        color_id = 0  # reset the color form the beginning
    
    #combine results
    
    #combined = combine_series(group_list, sigmathresh=sigmathresh)
    combined = combine_series_methodB(group_list, sigmathresh=sigmathresh)

    print("Final combination for marginalisation ", config.name_marg_spline)

    savefig = indiv_marg_dir + config.name_marg_spline + "_sigma_%2.2f_myplot.pdf" % sigmathresh
    
    #create plot
    stnd_plot.dmagplot(group_list,savefig,colors=colors,refgroup=combined)
    
    print("Saved group_list as ",str(marginalisation_dir +  config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl'),\
        "and combined result as ",str(marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl') )
    pkl.dump(group_list,
             open(marginalisation_dir  +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_groups.pkl',
                  'wb'))
    pkl.dump(combined,
             open(marginalisation_dir  +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % sigmathresh + '_combined.pkl',
                  'wb'))



# In[ ]:


if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Get the flux ratio from the observed lightcurves",
                               formatter_class=ap.RawTextHelpFormatter)
    help_lensname = "name of the lens to process"
    help_dataname = "name of the data set to process"
    parser.add_argument(dest='lensname', type=str,
                        metavar='lens_name', action='store',
                        help=help_lensname)
    parser.add_argument(dest='dataname', type=str,
                        metavar='dataname', action='store',
                        help=help_dataname)
    parser.add_argument('-v','--verbose',help="Verbosity",
                        dest="verbose", 
                        default=False,action="store_true")
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname

    present_program(sys.argv[0])
    sys.path.append("myconfig/")
    config_file = "myconfig_" + lensname + "_" + dataname

    config = importlib.import_module(config_file)
    # outdated
    #analysis_path = pth.Path(config.analysis_directory)
    #savefig_path = analysis_path/str("FluxRatio")
    #savefig_path.mkdir()
    #kw_FR = FR_kw(config,analysis_path,savefig_path)
    #plot_FR_kw(config,kw_FR,savefig_path)

    ### WIP ###
    combine_models_mag(config,sigmathresh=config.sigmathresh)
    
    success(sys.argv[0])

