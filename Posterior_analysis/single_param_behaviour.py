# 20-01-2023
# plot a single param behaviour over pso and mcmc

from Utils.tools import *
from Utils.get_res import *
from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as plt
from corner import quantile
from Plots.plotting_tools import averaged_plot
if __name__=="__main__":
    ############################
    present_program(sys.argv[0])
    ############################
    parser = ArgumentParser(description="Plot behaviour of the one single parameter over PSO and MCMC")
    parser.add_argument("-p","--parameter",type=str,dest="param", help="Which parameter to consider")
    parser.add_argument('SETTING_FILES',nargs="+",default=[],help="Setting file(s) to model")
    args          = parser.parse_args()
    param         = args.param
    setting_names = [s.replace(".py","") for s in get_setting_name(args.SETTING_FILES)]
    backup_path   = "backup_results"
    for setting_name in setting_names:
        setting       = get_setting_module(setting_name).setting()
        savemcmc_path = get_savemcmcpath(setting,backup_path)
        savefig_path  = get_savefigpath(setting, backup_path)
        param_mcmc    = get_mcmc_prm(setting,backup_path)
        try:
            _ = param_mcmc.index(param)
        except ValueError:
            raise RuntimeError(f"{param} is not within the MCMC parameters of {setting_name}:\n{param_mcmc}")
        pos_mcmc      = get_mcmc_smpl_for_prm(setting,param,backup_path=backup_path)
        lkl_mcmc      = get_mcmc_logL(setting_name=setting,backup_path=backup_path)
        pso_chain     = get_pso_chain(setting_name=setting_name,backup_path=backup_path)
        name_pso,steps_pso,param_pso = pso_chain
        try:
            i_prm = param_pso.index(param)
        except ValueError:
            raise RuntimeError(f"{param} is not within the PSO parameters of {setting_name}:\n{param_pso}")
        lkl,pos,vel = steps_pso
        pos_pso = np.transpose(pos)[i_prm].T.tolist()
        #vel_pso =  np.transpose(vel)[i_prm].T
        lkl_pso = lkl
        """
        combined_pos   = [*pos_pso,*pos_mcmc.tolist()]
        combined_lkl   = [*lkl_pso,*lkl_mcmc.tolist()]
        combined_steps = np.arange(0,len(combined_lkl))
        fig, ax = plt.subplots(2)
        fig.suptitle(f"Param {param}")
        ax[0].plot(combined_steps,combined_pos,color="r")
        ax[0].axvline(len(pos_pso),c="k")
        ax[0].set_xlabel("Steps (PSO+MCMC)")
        ax[0].set_ylabel("Position")
        ax[1].set_yscale("log")
        ax[1].plot(combined_steps,-np.array(combined_lkl),color="g")
        ax[1].axvline(len(pos_pso),c="k")
        ax[1].set_xlabel("Steps (PSO+MCMC)")
        ax[1].set_ylabel("-LogLikelihood")
        """
        fig, ax = plt.subplots(4,figsize=(18,12))
        fig.suptitle(f"Behaviour of Param {param}")
        steps_pso = np.arange(0,len(pos_pso))
        ax[0].plot(steps_pso,pos_pso,color="r") #= averaged_plot(ax[0],pos_pso,col="r",plot_scatter=False,num_average=20)
        pso_fin = pos_pso[-1]
        ax[0].axhline(pso_fin,label=f"final val:{np.round(pso_fin,3)}",c="k")
        ax[0].set_ylabel("Position")
        ax[0].legend()
        ax[1].set_yscale("log")
        ax[1].plot(steps_pso,-np.array(lkl_pso),color="g")
        ax[1].set_xlabel("Steps (PSO)")
        ax[1].set_ylabel("-LogLikelihood")
        steps_mcmc = np.arange(0,len(pos_mcmc))
        #ax[2].plot(steps_mcmc,pos_mcmc,color="b") # to change as in mcmc_behaviour_plot
        ax[2] = averaged_plot(ax[2],pos_mcmc,col="b",plot_scatter=True,num_average=200)
        mcmc_min,mcmc_med,mcmc_max = quantile(pos_mcmc,q=[0.16,0.5,0.84])
        #ax[2].fill_between(steps_mcmc,mcmc_min,mcmc_max,alpha=.5,color="g")
        ax[2].axhline(pso_fin,label=f"PSO final val:{np.round(pso_fin,3)}",c="k")
        err_mcmc = np.round(np.mean([mcmc_med-mcmc_min,mcmc_max-mcmc_med]),3)
        ax[2].axhline(mcmc_med,label=f"MCMC mean val:{np.round(mcmc_med,3)}+-{err_mcmc}",c="k",ls="--")
        ax[2].set_ylabel("Position")
        ax[2].legend()
        #ax[3].set_yscale("log")
        ax[3]=averaged_plot(ax[3],-np.array(lkl_mcmc),col="g",plot_scatter=False,num_average=200) #.plot(steps_mcmc,-np.array(lkl_mcmc),color="g")
        ax[3].set_xlabel("Steps (MCMC)")
        ax[3].set_ylabel("-LogLikelihood")
        image_path = create_path_from_list([savefig_path,f"{param}_bhv.png"])
        plt.tight_layout()
        plt.legend()
        plt.savefig(image_path)
        print(f"Saving image {image_path}")
    success(sys.argv[0])
