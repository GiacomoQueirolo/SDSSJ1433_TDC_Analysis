import pickle
import numpy as np
import matplotlib.pyplot as plt
import pathlib as pth
import sys,os
from corner import quantile
import importlib
import argparse as ap
import copy

from Utils.tools import *    
from stnd_handling_data import Error,Error_mag
from stnd_plot import print_res_w_err

plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['font.size']= 16

bins_n = 30

def plot_result(config,overplot=False,verbose=False):
    config = get_config(config)
    savefig_path=pth.Path(config.analysis_directory)
    if overplot:
        fgsize = 18
        fg,ax= plt.subplots(2,2,figsize=(fgsize,fgsize))
        fg_mag,ax_mag= plt.subplots(2,2,figsize=(fgsize,fgsize))
        colors = config.colors
        color_id = 0
        legends = []
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("analysis_kn"+str(knt))
        mkdir(savefig_path)        
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):
                mllist_name,mlfp = ml_config 
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                elif mltype_i=="splml":
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                        
                saveml_path = saveknt_path/config.get_savemlpath(mltype_i,ml_config)# saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str) 

                with open(str(saveml_path)+"/td.data", 'rb') as f:
                    timedelays=pickle.load(f)
                with open(str(saveml_path)+"/dmags.data", 'rb') as f:
                    magshifts=pickle.load(f)
                name  = str(config.combkw[a][mlt_i][mlc_i]).replace("spl1_","")
                distr = np.transpose(timedelays)
                distr_mag = np.transpose(magshifts)
                plot_result_single(distr,name,saveml_path)
                if verbose:
                    print("Done Resulting Distr "+name.replace("_"," "))
                plot_result_single(distr_mag,name,saveml_path,mag=True)                    
                if verbose:
                    print("Done Resulting Distr Mag "+name.replace("_"," "))
                if overplot:
                    ax,lgi   = overplot_result_single(distr,name,saveml_path,ax=ax,alpha=.5,color=colors[color_id])
                    ax_mag,_ = overplot_result_single(distr,name,saveml_path,ax=ax_mag,alpha=.5,color=colors[color_id])                    
                    color_id+=1
                    if color_id==len(colors):
                        color_id=0
                    legends.append(lgi)

    if overplot:

        for ax_del in (ax[0][1],ax_mag[0][1]):
            if len(legends)<20:
                ax_del.legend()   
            else:
                set_lgnd = ax_del.legend(handles=legends[:int(len(legends)/2.)],loc="upper left")
                ax_del.add_artist(set_lgnd)
                ax_del.legend(handles=legends[int(len(legends)/2.):],loc="upper center")
        
        
        overplot_name_dt  = f"{savefig_path}/Overplot_DtRes.pdf"
        overplot_name_mag = f"{savefig_path}/Overplot_MagRes.pdf"
        fg.savefig(overplot_name_dt)
        print(f"Saving {overplot_name_dt}") 
        fg_mag.savefig(overplot_name_mag)
        print(f"Saving {overplot_name_mag}") 
        

def overplot_result_single(distr,name,save_path,ax,alpha=.5,mag=False,color=None):
    # OverPlot the resulting time delay distribution
    if not mag:
        Dunit = "\Delta t"
        unit  = "d"
    else:
        Dunit = "\Delta\ mag"
        unit  = ""
    for i in range(len(distr)): #ABCD
        if i==0:           
            ax_i=ax[0][0] 
            let = "AB" 
        if i==1: 
            ax_i=ax[1][0] 
            let = "AC" 
        if i==2: 
            ax_i=ax[1][1]             
            let = "BC"            
        ddt = distr[i]
        arr,_,_ = ax_i.hist(ddt,bins=bins_n,alpha=alpha,color=color)
        ax_i.set_xlabel("$"+Dunit+"$ "+let+" ["+unit+"]")
        ax_i.set_yticklabels([])
        ax_i.set_yticks([])
    axdel=ax[0][1]
    lgd = axdel.plot(1,1,label=name.replace("_"," "),color=color)
    #axdel.legend()
    axdel.axis("off")
    return ax,lgd



def plot_result_single(distr,name,save_path,wD=False,mag=False):
    # Plot the resulting time delay distribution
    fg,ax= plt.subplots(2,2,figsize=(12,12))
    if not mag:
        Dunit = "\Delta t"
        unit  = "d"
    else:
        Dunit = "\Delta\ mag"
        unit  = ""
    plt.suptitle("Resulting $"+Dunit+"$ Distribution for "+name.replace("_"," "))
    for i in range(len(distr)): #ABCD
        if i==0:           
            ax_i=ax[0][0] 
            let = "AB" 
        if i==1: 
            ax_i=ax[1][0] 
            let = "AC" 
        if i==2: 
            ax_i=ax[1][1]             
            let = "BC"       
            if wD:
                let="AD"     
        ddt = distr[i]
        vl,vc,vr = quantile(ddt,q=[0.16,0.5,0.84])
        err_min = vc-vl
        err_max = vr-vc

        arr,_,_ = ax_i.hist(ddt,bins=bins_n,)
        label_res = "$"+Dunit+"\ "+str(np.round(vc,1))+"_{-"+str(np.round(err_min,1))+"}^{+"+str(np.round(err_max,1))+"} $"+unit 
        ax_i.axvline(vc,c="r",ls="--",label=label_res,linewidth=2.5) 
        max_plot = max(arr)*1.1
        ax_i.set_ylim(0,max_plot)
        ax_i.set_xlabel("$"+Dunit+"$ "+let+" ["+unit+"]")
        ax_i.fill_between(np.linspace(vl,vr) , 0, max_plot, color='grey', alpha=0.2)
        
        ax_i.legend()
    axdel=ax[0][1]
    axdel.plot(1,1,label=name.replace("_"," "),color="w")
    axdel.plot(1,1,label=r"Median",color="r",ls="--",linewidth=2.5)
    axdel.fill_between([1,1],1,1,color='grey', alpha=0.2, label="1-$\sigma$ region")
    axdel.legend()
    axdel.axis("off")
    if mag:
        name = "Mag_"+name
    plt.savefig(str(save_path)+"/ResDist_"+name+".png")
    plt.close()


def plot_err(config,verbose=False):    
    output_dir = config.lens_directory 
    opt = config.optset[0]
    base_lcs = config.get_lcs() 
    plot_path = config.figure_directory
    savefig_path=pth.Path(config.simu_directory)
    for knt_i, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        saveknt_path = savefig_path/str("simulation_kn"+str(knt))
        for mlt_i, mltype_i in enumerate(config.mltype): #ml type
            for mlc_i,ml_config in enumerate(config.ml_config[mlt_i]):  
                lcs = copy.deepcopy(base_lcs)
                mllist_name,mlfp = ml_config 
                """
                kwargs_ml = {"mltype":mltype_i,"mllist_name":mllist_name}
                if mltype_i=="polyml":
                    mlfp_str="_mlfp_"+str(mlfp)
                    kwargs_ml["mlfp"] = mlfp
                elif mltype_i=="splml":
                    if config.forcen:
                        mlfp_str="_nmlspl_"+str(mlfp)
                    else:
                        mlfp_str="_knstml_"+str(mlfp)
                    kwargs_ml["forcen"] = config.forcen
                    kwargs_ml["nmlspl"] = mlfp
                """                    
                combdir =  saveknt_path/config.get_savemlpath(mltype_i,ml_config)# saveknt_path/str("ml"+mltype_i[:-2]+mllist_name+mlfp_str)         

                #MOD_ROB
                if not check_success_analysis(get_analdir(combdir)):
                    print("Analysis from "+get_analdir(combdir)+" was not sufficiently precise. Ignored.")
                    continue

                sim_path  = str(combdir)+"/sims_"+config.simset_mock+"_opt_"+opt
        
                error = Error(sim_path)
                error_mag = Error_mag(sim_path)
                
                err_distr = error.get_distr() #not corrected for systematic
                err_mag_distr = error_mag.get_distr() #not corrected for systematic
                name = str(config.combkw[knt_i][mlt_i][mlc_i]).replace("spl1_","")
                plot_err_distr_single(err_distr,name,plot_path)
                if verbose:
                    print("Done Err Distr "+name.replace("_"," "))
                plot_err_distr_single(err_mag_distr,name,plot_path,mag=True)
                if verbose:
                    print("Done Err Distr Mag "+name.replace("_"," "))

def plot_err_distr_single(distr,name,plot_path,mag=False):
    fg,ax= plt.subplots(2,2,figsize=(12,12))
    plt.suptitle("Error Distribution of "+name.replace("_"," "))
    for i in range(len(distr)): 
        if i==0:           
            ax_i=ax[0][0] 
            let = "AB" 
        if i==1: 
            ax_i=ax[1][0] 
            let = "AC" 
        if i==2: 
            ax_i=ax[1][1]             
            let = "BC"            
        ddt = distr[i]
        vl,vc,vr = quantile(ddt,q=[0.16,0.5,0.84])
        err_min = vc-vl
        err_max = vr-vc

        sys_err = vc
        rnd_err = (err_min + err_max)/2.
        
        arr,_,_ = ax_i.hist(ddt,bins=bins_n,color="red")
        ax_i.axvline(vc,c="g",ls="--",linewidth=2.5)
        max_plot = max(arr)*1.1
        ax_i.set_ylim(0,max_plot)
        if not mag:
            unit  = "d"
            ax_i.set_xlabel("$\Delta t_{sim}-\Delta t_{meas}$ "+let+" ["+unit+"]")
        else:
            unit  = ""
            ax_i.set_xlabel("$\Delta mag_{sim}-\Delta mag_{meas}$ "+let+" ["+unit+"]")
        ax_i.fill_between(np.linspace(vl,vr) , 0, max_plot, color='grey', alpha=0.2)
        ax_i.axvline(0,c="k",ls="-.",linewidth=2.5)
        to_annotate = "      Error\n"
        ####
        ####
        sys_err_str = print_res_w_err(sys_err,sys_err).split("$\pm$")[0]
        rnd_err_str = print_res_w_err(rnd_err,rnd_err).split("$\pm$")[0]
        tot_err = np.sqrt(sys_err**2+ rnd_err**2)
        tot_err_str = print_res_w_err(tot_err,tot_err).split("$\pm$")[0]
        
        to_annotate += "Sys.: "+sys_err_str+" "+unit+""+"\n"+"Rnd.: "+rnd_err_str+" "+unit
        to_annotate += "\nTot.: "+tot_err_str+" "+unit
        #to_annotate += "Sys.: "+str(np.round(sys_err,2))+" "+unit+""+"\n"+"Rnd.: "+str(np.round(rnd_err,2))+" "+unit
        #to_annotate += "\nTot.: "+str(np.round(np.sqrt(sys_err**2+ rnd_err**2),2))+" "+unit
        ax_i.annotate(to_annotate,xy=(.65,.73),xycoords="axes fraction",bbox=dict(boxstyle='round', fc='w'))
    axdel=ax[0][1]
    axdel.plot(1,1,label=name.replace("_"," "),color="w")
    axdel.plot(1,1,label=r"Median",color="g",ls="--",linewidth=2.5)
    axdel.plot(1,1,label=r"Zero",color="k",ls="-.",linewidth=2.5)
    axdel.fill_between([1,1],1,1,color='grey', alpha=0.2, label="1-$\sigma$ region")
    axdel.legend()
    axdel.axis("off")
    if not mag:
        plt.savefig(str(plot_path)+"/ErrDist_"+name+".png")
    else:
        plt.savefig(str(plot_path)+"/ErrDistMag_"+name+".png")
    
    plt.close()

if __name__ == '__main__':
    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Plot distribution of time delay and Dmag analysis's result and/or error",
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
    parser.add_argument('-res','--result_distrib',help="Plot the time delay result distribution",
                        dest="res", 
                        default=False,action="store_true")
    parser.add_argument('-err','--error_distrib',help="Plot the time delay error distribution",
                        dest="err", 
                        default=False,action="store_true")
    args = parser.parse_args()
    lensname = args.lensname
    dataname = args.dataname
    verbose  = args.verbose
    res      = args.res
    err      = args.err
    
    name_prog = sys.argv[0]
    present_program(name_prog)
    if not res and not err:
        # if no input, it must be only the result distribution by default
        res = True
    
    sys.path.append("myconfig/")
    config_file = "myconfig_" + lensname + "_" + dataname
    config = importlib.import_module(config_file)    

    if res:
        plot_result(config,verbose=verbose)
    if err:
        plot_err(config,verbose=verbose)        
    success(name_prog)

