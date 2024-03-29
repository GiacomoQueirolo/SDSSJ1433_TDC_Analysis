## little functions to easily get the results from previous calculations
import json,pickle ,dill
from corner import quantile
from numpy import abs as npabs
from numpy import array

from Utils.tools import *

def load_whatever(name):
    try:        
        with open(name, 'r') as f:
            data = json.load(f)
    except:
        try:
            with open(name, 'rb') as f:
                data = pickle.load(f)
        except:
            try:
                with open(name,'rb') as f:
                    data = dill.load(f)
            except:
                with open(name, 'r') as f:
                    data = f.readlines()
                data = [data_l.replace(",\n","") for data_l in data]   
    return data

def get_kwres(setting_name,updated=False,backup_path="backup_results"):
    setting_name  = get_setting_name(setting_name)
    if not updated:
        savefig_path = get_savefigpath(setting_name,backup_path)
        name =  savefig_path+'/read_results.data'
    else:
        savemcmc_path = get_savemcmcpath(setting_name,backup_path)
        name =  savemcmc_path+'/read_results_updated.data'
    kwargs_result = load_whatever(name)
    return {"kwargs_results":kwargs_result}


def get_mcmc_prm(setting_name,backup_path="backup_results"):  
    setting_name  = get_setting_name(setting_name)
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
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
    
    return array(mcmc)
    
def get_mcmc_fermat(setting_name,backup_path="backup_results"):
    return get_mcmc_chain(setting_name,"ordered_fermat",backup_path)

"""
def get_mcmc_Df(setting_name,backup_path="backup_results",noD=True):
    # noD: ignore image D and return AB,AC and BC instead
    # return : mcmc_Df, shape: len_mcmc, 3
    mcmc_fermat = get_mcmc_fermat(setting_name,backup_path)
    # first MUST be A
    mcmc_DfT = np.transpose(mcmc_fermat)[1:]-np.transpose(mcmc_fermat)[0] 
    #mcmc_Df = mcmc_DfT.T.tolist() 
    if noD:
        mcmc_cp    = np.array(deepcopy(mcmc_DfT))
        # BC = C - B = (C-A)-(B-A) = AC - AB
        mcmc_BC    = mcmc_cp[1] - mcmc_cp[0]  
        mcmc_cp[2] = mcmc_BC
        mcmc_Df    = mcmc_cp.T.tolist()
    else:
        mcmc_Df    = mcmc_DfT.T.tolist()
    return mcmc_Df    
"""

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
    

def get_pso_chain(setting_name,backup_path="backup_results",use_backup=False): 
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    try:
        if use_backup:
            raise FileNotFoundError("use_backup=True")
        pso_file_name = save_json_name(setting_name,savemcmc_path,"pso")
        pso_chain     = load_whatever(pso_file_name)
    except FileNotFoundError:
        pso_chain     = load_whatever(f"{savemcmc_path}/psobackup.json")
    return pso_chain
    
def get_mcmc(setting_name,backup_path="backup_results"):
    mcmc_smpl = get_mcmc_smpl(setting_name,backup_path)
    mcmc_prm  = get_mcmc_prm(setting_name,backup_path)
    mcmc_logL = get_mcmc_logL(setting_name,backup_path)
    return {"mcmc_smpl":mcmc_smpl,"mcmc_prm":mcmc_prm, "mcmc_logL":mcmc_logL}

"""    
def get_results(setting_name,backup_path="backup_results"):
    mcmc_res = get_mcmc(setting_name,backup_path)
    kw_res   = get_kwres(setting_name,backup_path=backup_path)
    return {**mcmc_res,**kw_res}
"""


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



def get_sigma_kw(setting,mcmc_chain=None,print_res=None,save=True,backup_path="./backup_results/"):
    sett = get_setting_module(setting,1)
    if mcmc_chain is None:
        _mcmc_chain = get_mcmc(sett).values()
        mcmc_chain  = ["MCMC",*_mcmc_chain]
    kwargs_sigma_upper={}
    kwargs_sigma_lower={}
    n_ra,n_dec = 0,0
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = mcmc_chain 
    if print_res is None:
        print_res = open(get_savefigpath(sett,backup_path)+"/results.txt","a")
        
    for i in range(len(param_mcmc)):
        val_min, val, val_max = quantile(samples_mcmc[:,i],q=[0.16, 0.5, 0.84])
        sig_min  = npabs(val_min-val)
        sig_max  = val_max - val
        if print_res is not False:
            str_val  = "{:.2e}".format(val)
            str_min  = "{:.2e}".format(sig_min)
            str_max  = "{:.2e}".format(sig_min)
            if str_max==str_min:
                print_res.write(param_mcmc[i]+" "+ str_val+" +-  "+str_min+"\n")
            else:
                print_res.write(param_mcmc[i]+" "+str_val+" - "+str_min+" + "+str_max+"\n")
        if param_mcmc[i]!="ra_image" and param_mcmc[i]!="dec_image":
            kwargs_sigma_lower[param_mcmc[i]]=sig_min
            kwargs_sigma_upper[param_mcmc[i]]=sig_max

        elif param_mcmc[i]=="ra_image":
            kwargs_sigma_lower["ra_image_"+str(n_ra)]=sig_min
            kwargs_sigma_upper["ra_image_"+str(n_ra)]=sig_max            
            n_ra+=1
        else:
            kwargs_sigma_lower["dec_image_"+str(n_dec)]=sig_min
            kwargs_sigma_upper["dec_image_"+str(n_dec)]=sig_max            
            n_dec+=1
    if print_res is not False:
        print_res.write("\n#################################\n")
    if save:
        for name,res in zip(["read_sigma_low","read_sigma_up"],[kwargs_sigma_lower,kwargs_sigma_upper]):
            save_name = get_savefigpath(sett,backup_path)+"/"+name+".data"
            print("Saving "+save_name)
            pickle.dump(res,open(save_name,"wb"))
    else:
        return {"read_sigma_low":kwargs_sigma_lower,"read_sigma_up":kwargs_sigma_upper}
    

def get_combined_Df(combined_sett,main_dir="./"):
    combined_setting = get_combined_setting_module(combined_sett,main_dir=main_dir)
    dir_path         = combined_setting.get_respath()
    KDE              = combined_setting.KDE
    
    Combined_PDF = load_whatever(f"{dir_path}/Combined_PDF{['_KDE' if KDE else ''][0]}.pkl")
    
        
    if KDE:
        Positions_KDE = load_whatever(f"{dir_path}/Combined_PDF_KDE_positions.pkl")
    else:
        Combined_bins = load_whatever(f"{dir_path}/Combined_PDF_bins.pkl")

    if KDE:
        return combined_setting, np.array(Combined_PDF),np.array(Positions_KDE)
    else:
        return combined_setting, np.array(Combined_PDF),np.array(Combined_bins)
    
    

def get_combined_Rmag(combined_sett,main_dir="./"):
    combined_setting = get_combined_setting_module(combined_sett,main_dir=main_dir)
    dir_path_mag     = combined_setting.get_respath()+"/Mag/"
    KDE              = combined_setting.KDE

    Combined_mag_PDF = load_whatever(f"{dir_path_mag}/Combined_mag_PDF{['_KDE' if KDE else ''][0]}.pkl")
        
    if KDE:
        Positions_mag_KDE = load_whatever(f"{dir_path_mag}/Combined_mag_PDF_KDE_positions.pkl")
    else:
        Combined_mag_bins = load_whatever(f"{dir_path_mag}/Combined_mag_PDF_bins.pkl")
        
    if KDE:
        return combined_setting, np.array(Combined_mag_PDF),np.array(Positions_mag_KDE)
    else:
        return combined_setting, np.array(Combined_mag_PDF),np.array(Combined_mag_bins)
    
