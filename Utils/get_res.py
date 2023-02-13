## little functions to easily get the results from previous calculations
from Utils.tools import *
import json,pickle 
from corner import quantile
from numpy import abs as npabs

from lenstronomy.Sampling.parameters import Param
from Data.input_data import get_kwargs_model,get_kwargs_constraints

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
    savemcmc_patdef get_prm_list(setting):
    return get_Param(setting).num_param()[1]
h = get_savemcmcpath(setting_name,backup_path)
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
    

def get_pso_chain(setting_name,backup_path="backup_results"): 
    setting_name  = get_setting_name(setting_name)
    savemcmc_path = get_savemcmcpath(setting_name,backup_path)
    file_name     = setting_name.replace("settings","pso").replace(".py","")+".json"
    pso_file_name = create_path_from_list([savemcmc_path,file_name])
    pso_chain     = load_whatever(pso_file_name)
    return pso_chain
    def get_prm_list(setting):
    return get_Param(setting).num_param()[1]

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


def conv_mcmc_i_to_kwargs(setting,mcmc_i):
    param_class   = get_Param(setting)
    kwargs_result = param_class.args2kwargs(mcmc_i, bijective=True)
    return kwargs_result

def get_Param(setting):
    setting = get_setting_module(setting,1)
    kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens = setting.lens_params
    kwargs_ps_init, kwargs_ps_sigma, kwargs_fixed_ps, kwargs_lower_ps, kwargs_upper_ps = setting.ps_params
    kwargs_lens_light_init, kwargs_lens_light_sigma, kwargs_fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light = setting.lens_light_params
    if not setting.WS:
        kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source = setting.source_params
    else:
        kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source = None,None,None,None,None
    
    kwargs_model       = get_kwargs_model(setting)
    kwargs_constraints = get_kwargs_constraints(setting)
    param_class = Param(kwargs_model, kwargs_fixed_lens=kwargs_fixed_lens, kwargs_fixed_source=kwargs_fixed_source,
                        kwargs_fixed_lens_light=kwargs_fixed_lens_light, kwargs_fixed_ps=kwargs_fixed_ps, kwargs_fixed_special=None,
                        kwargs_fixed_extinction=None, 
                        kwargs_lower_lens=kwargs_lower_lens, kwargs_lower_source=kwargs_lower_source, 
                        kwargs_lower_lens_light=kwargs_lower_lens_light, kwargs_lower_ps=kwargs_lower_ps,
                        kwargs_lower_special=None, kwargs_lower_extinction=None,
                        kwargs_upper_lens=kwargs_upper_lens, kwargs_upper_source=kwargs_upper_source, 
                        kwargs_upper_lens_light=kwargs_upper_lens_light, kwargs_upper_ps=kwargs_upper_ps,
                        kwargs_upper_special=None, kwargs_upper_extinction=None,
                        kwargs_lens_init=None, **kwargs_constraints)
    return param_class

def get_prm_list(setting,backup_path="./backup_path"):
    list_prm = get_Param(setting).num_param()[1]
    # validity check:
    try:
        list_prm_mcmc=get_mcmc_prm(setting,backup_path=backup_path)
        if list_prm!=list_prm_mcmc:
            raise RuntimeError("The parameter have changed since the MCMC run!")
    except FileNotFoundError:
        print("warning: I couldn't double check that the parameter didn't change since the MCMC run")
    return list_prm


def get_sigma_kw(setting,mcmc_chain=None,print_res=None,save=True):
    sett = get_setting_module(setting,1)
    if mcmc_chain is None:
        _mcmc_chain = get_mcmc(sett).values()
        mcmc_chain  = ["MCMC",*_mcmc_chain]
    kwargs_sigma_upper={}
    kwargs_sigma_lower={}
    n_ra,n_dec = 0,0
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = mcmc_chain 
    if print_res is None:
        print_res = open(get_savefigpath(sett)+"/results.txt","a")
        
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
            kwargs_sigma_lower["ra_image_"+plottingstr(n_ra)]=sig_min
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
            save_name = get_savefigpath(sett)+"/"+name+".data"
            print("Saving "+save_name)
            pickle.dump(res,open(save_name,"wb"))
    else:
        return {"read_sigma_low":kwargs_sigma_lower,"read_sigma_up":kwargs_sigma_upper}
    
