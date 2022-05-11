## little functions to easily get the results from previous calculations
from tools import *
import json,pickle

def load_whatever(name):
    try:        
        with open(name, 'r') as f:
            data = json.load(f)
    except:
        try:
            with open(name, 'rb') as f:
                data = pickle.load(f)
        except:
            with open(name, 'r') as f:
                data = f.readlines()
            data = [data_l.replace(",\n","") for data_l in data]   
    return data

def get_kwres(setting_name,updated=False,backup_path="backup_results"):
    setting_name  = get_setting_name(setting_name)
    savefig_path = get_savefigpath(setting_name,backup_path)
    if not updated:
        name =  savefig_path+'/read_results.data'
    else:
        name =  savefig_path+'/read_results_updated.data'
    kwargs_result = load_whatever(name)
    return {"kwargs_results":kwargs_result}


def get_mcmc_prm(setting_name,backup_path="backup_results"):  
    setting_name  = get_setting_name(setting_name)
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    setting_name  = get_setting_name(setting_name)
    file_name = setting_name.replace("settings","mcmc_prm").replace(".py","")+".dat"
    mcmc_prm_file_name = savemcmc_path+file_name
    mcmc_prm = load_whatever(mcmc_prm_file_name)
    return mcmc_prm
    
def get_mcmc_chain(setting_name,mcmc_name,backup_path="backup_results"):
    setting_name  = get_setting_name(setting_name)
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    not_in =True
    for nm in ["smpl","logL","ordered_fermat","mag_rt"]:
        if mcmc_name ==nm:
            not_in=False
            break
    if not_in:
        raise RuntimeWarning("Give correct mcmc name, not "+mcmc_name)
    file_name      = setting_name.replace("settings","mcmc_"+mcmc_name).replace(".py","")+".json"
    mcmc_file_name = savemcmc_path+file_name
    mcmc = load_whatever(mcmc_file_name)
    from numpy import array
    return array(mcmc)
    
def get_mcmc_fermat(setting_name,backup_path="backup_results"):
    return get_mcmc_chain(setting_name,"ordered_fermat",backup_path)

def get_mcmc_mag(setting_name,backup_path="backup_results"):
    return get_mcmc_chain(setting_name,"mag_rt",backup_path)

def get_mcmc_logL(setting_name,backup_path="backup_results"):
    return get_mcmc_chain(setting_name,"logL",backup_path)
    
def get_mcmc_smpl(setting_name,backup_path="backup_results"):  
    return get_mcmc_chain(setting_name,"smpl",backup_path)  

def get_mcmc_smpl_for_prm(setting_name,prm_name,backup_path="backup_results"):  
    smpl = get_mcmc_smpl(setting_name,backup_path)  
    all_prm   = get_mcmc_prm(setting_name,backup_path=backup_path)
    try:
        index_prm = all_prm.index(prm_name)
    except ValueError:
        print(prm_name," not found, available params are:\n",all_prm)
    return smpl[:,index_prm]
    
    
    
    
def get_mcmc(setting_name,backup_path="backup_results"):
    mcmc_smpl = get_mcmc_smpl(setting_name,backup_path)
    mcmc_prm  = get_mcmc_prm(setting_name,backup_path)
    mcmc_logL = get_mcmc_logL(setting_name,backup_path)
    return {"mcmc_smpl":mcmc_smpl,"mcmc_prm":mcmc_prm, "mcmc_logL":mcmc_logL}

    
def get_results(setting_name,backup_path="backup_results"):
    mcmc_res = get_mcmc(setting_name,backup_path)
    kw_res   = get_kwres(setting_name,backup_path)
    return {**mcmc_res,**kw_res}


def check_logL(setting_name,backup_path="backup_results",max_diff=5,index_logl=2000): 
    logL = get_mcmc_logL(setting_name,backup_path)
    accept = True 
    logL_0 = logL[-index_logl] 
    for logL_i in logL[-index_logl:]: 
        diff = abs(logL_i-logL_0)
        if diff<max_diff: 
            accepted=True 
        else: 
            accepted=False
            print("LogLikelihood difference:",diff) 
            break 
    return accepted 


