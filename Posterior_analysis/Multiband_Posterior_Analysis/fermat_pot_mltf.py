#!/usr/bin/env python
# coding: utf-8

#  I want to correct the images position, which now is fixed, to be set to the one of the corresponding
#  mcmc step

import sys
import corner
import argparse
import numpy as np
import multiprocess

from Utils.order_images import image_order_mltf
#from Utils.get_res import *
from Utils.Multiband_Utils.tools_multifilter import *
from Utils.Multiband_Utils.get_res_multifilter import *
from Data.Multiband_Data.Param_mltf import get_Param_mltf
from Data.Multiband_Data.input_data_mltf import init_lens_model
from Posterior_analysis.fermat_pot_analysis import labels_Fermat,labels_Df,get_Df_from_frm,_get_fermat,gen_mcmc_Df


  

def gen_mcmc_fermat_mltf(mcmc,multifilter_setting,lensModel=None,param_class=None,verbose=False,backup_path="backup_path_mltf"):
    # mcmc_prior shape = param, n_points
    # setting   = setting module
    if lensModel is None:
        lensModel = init_lens_model(multifilter_setting)
    mcmc_fermat = []
    if param_class is None:
        param_class = get_Param_mltf(multifilter_setting)
    def getferm(mcmc_i):
        return _get_fermat(mcmc_i,param_class,lensModel)
        
    with multiprocess.Pool() as pool:
        mcmc_fermat = pool.map(getferm,mcmc) #[_get_fermat(mcmc_i,param_class,lensModel) for mcmc_i in mcmc]
    #I want to obtain the correct image order
    ########################################
    # the order should be independent on the setting
    #new_order = get_new_image_order(multifilter_setting.settings[0],mcmc,verbose=verbose,backup_path=backup_path)
    new_order   = image_order_mltf(multifilter_setting,mcmc,verbose=verbose,backup_path=backup_path)
    
    tmp_mcmc    = np.array(mcmc_fermat).transpose()
    mcmc_fermat = np.transpose([tmp_mcmc[i] for i in new_order])
    return mcmc_fermat.tolist() # shape: (steps, f(i) )



def get_mcmc_Df_mltf(multifilter_setting ,noD=True):
    # noD: ignore image D and return AB,AC and BC instead
    # return : mcmc_Df, shape: len_mcmc, 3
    mcmc_fermat = get_mcmc_fermat_mltf(multifilter_setting)
    mcmc_Df     = get_Df_from_frm(np.transpose(mcmc_fermat),BC=noD).T.tolist()
    return mcmc_Df 



@check_mltf_setting
def save_Df_mltf(multifilter_setting,no_plot=False):
    print_setting(multifilter_setting.name)
    ############################################################
    #This should be the same for all settings
    if not multifilter_setting.CP:
        lens_model_list = ['SIE']
    else:
        print("WARNING: Considering the PEMD profile for the main lens")
        lens_model_list = ['PEMD']
    lens_model_list= [*lens_model_list,'SIS','SHEAR_GAMMA_PSI']

    ############################################

    #MCMC sample
    samples_mcmc = get_mcmc_smpl_mltf(multifilter_setting)
    mcmc_fermat  = gen_mcmc_fermat_mltf(samples_mcmc,multifilter_setting)
    mcmc_Df      = gen_mcmc_Df(samples_mcmc,multifilter_setting,mcmc_fermat=mcmc_fermat)

    #Save the mcmc in a file, NOTE: they are ordered A,B,C,D
    mcmc_file_name = multifilter_setting.get_savejson_name("mcmc_ordered_fermat")
    multifilter_setting.savejson_data(mcmc_fermat,mcmc_file_name)

    if not no_plot:
        mcmc_fermat = np.array(mcmc_fermat[cut_mcmc:])
        mcmc_Df     = np.array(mcmc_Df[cut_mcmc:])

        plot = corner.corner(mcmc_fermat, labels=labels_Fermat, show_titles=True)
        plot.savefig(f"{multifilter_setting.savefig_path}/Single_fermat_potential.png")

        plot = corner.corner(mcmc_Df, labels=labels_Df, show_titles=True)
        plot.savefig(f"{multifilter_setting.savefig_path}/Single_Df.png")




if __name__=="__main__":
    #############################
    present_program(sys.argv[0])
    #############################
    parser = argparse.ArgumentParser(description="Produces the stacked MCMC results for the fermat potential at the position of the images - MULTIFILTER VERSION")
    parser.add_argument("-c", "--cut_mcmc", type=int, dest="cut_mcmc", default=0,
                        help="cut the first <c> steps of the mcmc to ignore them")
    parser.add_argument("-NP", "--no_plot",dest="no_plot", default=False,action="store_true",
                        help="Ignore the corner plots")
    parser.add_argument('MULTIFILTER_SETTING',nargs="1",default=[],help="Multifilter setting file to consider")
    
    args = parser.parse_args()
    multifilter_sett =  get_multifilter_setting_module(args.MULTIFILTER_SETTING)
    cut_mcmc = int(args.cut_mcmc)
    no_plot  = args.no_plot 
    save_Df_mltf(multifilter_sett,no_plot)
    success(sys.argv[0])
