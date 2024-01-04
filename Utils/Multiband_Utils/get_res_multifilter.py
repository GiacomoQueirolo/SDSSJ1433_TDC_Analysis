## little functions to easily get the results from previous calculations
import pickle 
from copy import deepcopy
from corner import quantile
from numpy import array
from numpy import abs as npabs

from Utils.get_res import *
from Utils.Multiband_Utils.tools_multifilter import *



@check_mltf_setting
def get_kwres_mltf(multifilter_setting,updated=False):
    if not updated:
        savefig_path = multifilter_setting.savefig_path
        name =  savefig_path+'/read_results.data'
    else:
        savemcmc_path = multifilter_setting.savemcmc_path
        name =  savemcmc_path+'/read_results_updated.data'
    kwargs_result = load_whatever(name)
    return {"kwargs_results":kwargs_result}

def get_mcmc_prm_mltf(multifilter_setting):
    return get_mcmc_chain_mltf(multifilter_setting,"prm")
    
@check_mltf_setting
def get_mcmc_chain_mltf(multifilter_setting,mcmc_name):
    not_in =True
    for nm in ["smpl","logL","ordered_fermat","mag_rt","prm"]:
        if mcmc_name == nm:
            mcmc_name = "mcmc_"+mcmc_name
            not_in=False
            break
    if not_in:
        raise RuntimeWarning("Give correct mcmc name, not "+mcmc_name)
    mcmc_file_name = multifilter_setting.get_savejson_path(mcmc_name)
    if "prm" in mcmc_name:
        mcmc_file_name = mcmc_file_name.replace(".json",".dat")
    mcmc = load_whatever(mcmc_file_name)
    return array(mcmc)
    
def get_mcmc_fermat_mltf(multifilter_setting):
    return get_mcmc_chain_mltf(multifilter_setting,"ordered_fermat")

def get_mcmc_Df_mltf(multifilter_setting,noD=True):
    # noD: ignore image D and return AB,AC and BC instead
    # return : mcmc_Df, shape: len_mcmc, 3
    mcmc_fermat = get_mcmc_fermat_mltf(multifilter_setting)
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


def get_mcmc_mag_mltf(multifilter_setting):
    return get_mcmc_chain_mltf(multifilter_setting,"mag_rt")

def get_mcmc_logL_mltf(multifilter_setting):
    return get_mcmc_chain_mltf(multifilter_setting,"logL")
    
def get_mcmc_smpl_mltf(multifilter_setting):  
    return get_mcmc_chain_mltf(multifilter_setting,"smpl")  

def get_mcmc_smpl_for_prm_mltf(multifilter_setting,prm_name):  
    smpl    = get_mcmc_smpl_mltf(multifilter_setting)  
    all_prm = get_mcmc_prm_mltf(multifilter_setting)
    try:
        index_prm = all_prm.index(prm_name)
    except ValueError:
        print(prm_name," not found, available params are:\n",all_prm)
    return smpl[:,index_prm]
    
@check_mltf_setting
def get_pso_chain_mltf(multifilter_setting): 
    pso_file_name = multifilter_setting.get_savejson_path("pso")
    pso_chain     = load_whatever(pso_file_name)
    return pso_chain
    
def get_mcmc_mltf(multifilter_setting):
    mcmc_smpl = get_mcmc_smpl_mltf(multifilter_setting)
    mcmc_prm  = get_mcmc_prm_mltf(multifilter_setting)
    mcmc_logL = get_mcmc_logL_mltf(multifilter_setting)
    return {"mcmc_smpl":mcmc_smpl,"mcmc_prm":mcmc_prm, "mcmc_logL":mcmc_logL}


def check_logL_mltf(multifilter_setting,max_diff=5,index_logl=2000): 
    logL = get_mcmc_logL_mltf(multifilter_setting)
    accepted = True 
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



def get_sigma_kw_mltf(multifilter_setting,mcmc_chain=None,print_res=None,save=True):
    mltf_sett = get_multifilter_setting_module(multifilter_setting)
    if mcmc_chain is None:
        _mcmc_chain = get_mcmc_mltf(mltf_sett).values()
        mcmc_chain  = ["MCMC",*_mcmc_chain]
    kwargs_sigma_upper={}
    kwargs_sigma_lower={}
    n_ra,n_dec = 0,0
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = mcmc_chain 
    if print_res is None:
        print_res = open(f"{mltf_sett.savefig_path}/results.txt","a")
        
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
            save_name = f"{mltf_sett.savefig_path}/{name}.data"
            print("Saving "+save_name)
            pickle.dump(res,open(save_name,"wb"))
    else:
        return {"read_sigma_low":kwargs_sigma_lower,"read_sigma_up":kwargs_sigma_upper}
    
