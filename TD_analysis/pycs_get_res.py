# Similar to get_res.py used in lenstronomy
# should simplify the retrieval of data

import importlib
import sys
import pickle
from TD_analysis.stnd_handling_data import Error,Error_mag
from Utils.tools import get_config,check_success_analysis,get_analdir
analysis_data_name = ["in_time_shift","in_mag_shift","chi_red","resid","splines","td","kwargs_mc"]#tot_mag_shift is only on polyml

        
def get_single_analysis(path_analysis,knt,mltype,mlconfig):
    mlLC,mlDof = mlconfig
    if mltype=="polyml":
        mlstr ="poly"+str(mlLC)+"_mlfp_"+str(mlDof)
    elif mltype=="splml":
        mlstr ="spl"+str(mlLC)+"_nmlspl_"+str(mlDof)
    path = str(path_analysis)+"/analysis_kn"+str(knt)+"/ml"+mlstr
    kw_analysis = {"knt":knt,"mltype":mltype,"mlconfig":mlconfig}
    for i in range(len(analysis_data_name)):
        data = pickle.load(open(path+"/"+analysis_data_name[i]+".data","rb"))
        kw_analysis[analysis_data_name[i]] = data
    #if mltype=="polyml":
    #    kw_analysis["tot_mag_shift"] = pickle.load(open(path+"/tot_mag_shift.data","rb"))
    return kw_analysis
    
    
    
def get_analysis(_config):
    config = get_config(_config)
    all_analysis  = []
    path_analysis = config.lens_directory+"/Analysis/"
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        for mlt_i, mltype in enumerate(config.mltype): #ml type
            for mlc_i,mlconfig in enumerate(config.ml_config[mlt_i]):
                kw_res = get_single_analysis(path_analysis,knt,mltype,mlconfig)
                all_analysis.append(kw_res)
    return all_analysis
    
    
def get_single_sim(path_sim,knt,mltype,mlconfig,simset_mock,opt):
    mlLC,mlDof = mlconfig
    if mltype=="polyml":
        mlstr ="poly"+str(mlLC)+"_mlfp_"+str(mlDof)
    elif mltype=="splml":
        mlstr ="spl"+str(mlLC)+"_nmlspl_"+str(mlDof)
    path = str(path_sim)+"/simulation_kn"+str(knt)+"/ml"+mlstr
    sim_path  = str(path)+"/sims_"+simset_mock+"_opt_"+opt+"/"
    kw_sim = {"knt":knt,"mltype":mltype,"mlconfig":mlconfig,"simset_mock":simset_mock,"opt":opt}
    if not check_success_analysis(get_analdir(path)):
        print("Analysis from "+get_analdir(path)+" was not sufficiently precise. Ignored .") 
        return 0                    
    kw_sim["error"]     = Error(sim_path)
    kw_sim["error_mag"] = Error_mag(sim_path) 
    return kw_sim

def get_error(_config):
    config = get_config(_config)
    all_errors  = []
    path_sim = config.lens_directory+"/Simulation/"
    opt = config.optset[0]
    simset_mock = config.simset_mock
    for a, knt in enumerate(config.knotstep_marg): #knotstep of the intr.spline
        for mlt_i, mltype in enumerate(config.mltype): #ml type
            for mlc_i,mlconfig in enumerate(config.ml_config[mlt_i]):
                kw_sim = get_single_sim(path_sim,knt,mltype,mlconfig,simset_mock,opt)
                if kw_sim==0:
                    continue
                all_errors.append(kw_sim)
    return all_errors

def get_combined_res(_config,main_dir_path="."):
    config = get_config(_config,config_path=main_dir_path+"/myconfig/")
    marginalisation_dir = config.combined_directory + config.name_marg_spline + '/'
    combined = pickle.load(open(main_dir_path+"/"+marginalisation_dir + config.name_marg_spline + "_sigma_%2.2f" % config.sigmathresh + '_combined.pkl','rb'))
    return combined


def get_combined_mag(_config,main_dir_path="."):
    config = get_config(_config,config_path=main_dir_path+"/myconfig/")
    marginalisation_dir = config.combined_directory + config.name_marg_spline + '/'
    combined_mag = pickle.load(open(main_dir_path+"/"+marginalisation_dir +"mag_"+ config.name_marg_spline + "_sigma_%2.2f" % config.sigmathresh + '_combined.pkl','rb'))
    return combined_mag
    
