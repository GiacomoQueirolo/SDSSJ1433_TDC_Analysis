from lenstronomy.Sampling.parameters import Param

from Utils.tools import * 
from Utils.get_res import get_mcmc_prm
from Utils.Multiband_Utils.tools_multifilter import * 
from Data.Multiband_Data.input_data_mltf import get_kwargs_model_mltf, get_kwargs_constraints_mltf

@check_mltf_setting
def _get_Param_mltf(multifilter_sett):
    kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens = multifilter_sett.lens_params
    kwargs_ps_init, kwargs_ps_sigma, kwargs_fixed_ps, kwargs_lower_ps, kwargs_upper_ps = multifilter_sett.ps_params
    kwargs_lens_light_init, kwargs_lens_light_sigma, kwargs_fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light = multifilter_sett.lens_light_params
    if not multifilter_sett.allWS:
        kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source = multifilter_sett.source_params
    else:
        kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source = None,None,None,None,None
    
    kwargs_model       = get_kwargs_model_mltf(multifilter_sett)
    kwargs_constraints = get_kwargs_constraints_mltf(multifilter_sett,kwargs_model=kwargs_model)
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

def get_Param_mltf(multifilter_sett):#,save=True):
    # problem: can't pickle objects
    """try:
        with open(get_savefigpath(setting)+"/Prm_class.pkl","rb") as f:
            param_class = pickle.load(f)
    except FileNotFoundError:
        param_class = _get_Param(setting)
        if save:
            with open(get_savefigpath(setting)+"/Prm_class.pkl","wb") as f:
                pickle.dump(param_class,f)
    """
    return _get_Param_mltf(multifilter_sett)

def get_prm_list_mltf(multifilter_sett):
    param_class = get_Param_mltf(multifilter_sett)
    list_prm    = param_class.num_param()[1]
    # validity check:
    """
    try:
        list_prm_mcmc=get_mcmc_prm(setting,backup_path=backup_path)
        if list_prm!=list_prm_mcmc:
            raise RuntimeError("The parameter have changed since the MCMC run!")
    except FileNotFoundError:
        print("warning: I couldn't double check that the parameter didn't change since the MCMC run")
    """
    # for now ignored
    return list_prm

def conv_mcmc_i_to_kwargs_mltf(multifilter_sett,mcmc_i,Param_class=None):
    if Param_class is None:
        Param_class   = get_Param_mltf(multifilter_sett)
    kwargs_result = Param_class.args2kwargs(mcmc_i, bijective=True)
    return kwargs_result

@check_mltf_setting
def get_sigma_mltf(multifilter_setting):
    keywords = "lens","source","lens_light","ps","special","extinction"
    sigma_dic = {}
    for kw in keywords:
        not_present = []
        if kw == "special":
            not_present = {}
        sigma_dic[f"kwargs_{kw}"] = getattr(multifilter_setting,f"{kw}_params",[None,not_present])[1]
    return sigma_dic

# to implement further
"""
def count_images(params):
    n_im_ra  = 0
    n_im_dec = 0
    for p in params:
        if "ra_image" in p:
            n_im_ra+=1
        elif "dec_image" in p:
            n_im_dec+=1
    if n_im_dec!=n_im_ra:
        raise RuntimeError("Number of coordinates for images in parameters not matching")
    return n_im_ra

def index_common_params(sett_list,param_list=None,no_light=True):
    if param_list is None:
        param_list=[get_prm_list(s) for s in sett_list]
    common_params = set(param_list[0])
    for prm in param_list[1:]:
        common_params.intersection_update(prm)
    common_params=list(common_params)
    if no_light:
        # ignore light params
        common_params = [c for c in common_params if not "light" in c]
    # note: ra_image and dec_image are only once there
    indexes = [flatten([[ip for cp in common_params if pi==cp] for ip,pi in enumerate(p)]) for p in param_list]
    return indexes

    

"""